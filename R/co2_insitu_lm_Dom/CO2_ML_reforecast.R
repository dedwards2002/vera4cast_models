# ==============================================================================
# MACHINE LEARNING REFORECAST ANALYSIS: Random Forest CO2 vs Baselines
# ==============================================================================

# ------ Load packages -----
library(tidyverse)
library(lubridate)
library(vera4castHelpers)
library(tidymodels)
library(score4cast)
library(duckdb)
library(DBI)

# --- Configuration ---
my_model_id <- 'Dom_RF_Reforecast'
focal_sites <- c("fcre")
target_variable <- "CO2_umolL_sample"
covariate_variables <- c("Chla_ugL_mean", "CO2flux_umolm2s_mean") 
met_variables <- c("air_temperature")

# We will run 1 forecast per month for the year 2023 to save computation time
# (You can change this to "week" but Random Forest training takes longer than LM)
reforecast_dates <- seq(ymd('2023-01-01'), ymd('2023-12-31'), by = "month")
horizon_days <- 30 
n_members <- 31

#------- 1. Read Target Data --------
message("--- Downloading VERA Targets ---")
targets_url <- "https://amnh1.osn.mghpcc.org/bio230121-bucket01/vera4cast/targets/project_id=vera4cast/duration=P1D/daily-insitu-targets.csv.gz"
targets <- read_csv(targets_url, show_col_types = FALSE) |>
  filter(site_id %in% focal_sites, variable %in% c(target_variable, covariate_variables)) |>
  mutate(datetime = as_date(datetime))

# Extract the CO2 targets
targets_co2 <- targets |> 
  filter(variable == target_variable, depth_m == 0.1) |> 
  group_by(datetime, site_id) |> 
  summarize(co2 = mean(observation, na.rm = TRUE), .groups = "drop")

# Extract the covariates
covariates_historical <- targets |> 
  filter(variable %in% covariate_variables) |> 
  select(datetime, site_id, variable, observation) |> 
  group_by(datetime, site_id, variable) |> 
  summarize(observation = mean(observation, na.rm = TRUE), .groups = "drop") |> 
  pivot_wider(names_from = variable, values_from = observation)

# ------ 2. Download Historical Weather  ------
message("--- Downloading Historical NOAA Stage 3 ---")
weather_past <- vera4castHelpers::noaa_stage3() |> 
  dplyr::filter(site_id %in% focal_sites, datetime >= ymd('2018-01-01'), variable %in% met_variables) |> 
  dplyr::collect()

weather_past_daily <- weather_past |> 
  mutate(datetime = as_date(datetime)) |> group_by(datetime, site_id, variable) |> 
  summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
  mutate(prediction = prediction - 273.15) |> pivot_wider(names_from = variable, values_from = prediction)

# Build the complete historical driver dataset
model_drivers <- weather_past_daily |> 
  left_join(covariates_historical, by = c("datetime", "site_id")) |> 
  distinct(datetime, site_id, .keep_all = TRUE) |> 
  group_by(site_id) |> arrange(datetime) |> 
  tidyr::fill(all_of(covariate_variables), .direction = "down") |> ungroup()

# Combine targets and drivers for the training set
targets_rf <- targets_co2 |> 
  left_join(model_drivers, by = c("datetime", "site_id")) |> 
  arrange(site_id, datetime) |> group_by(site_id) |> 
  mutate(co2_yday = lag(co2, 1)) |> ungroup() |> 
  filter(complete.cases(co2, co2_yday, air_temperature, Chla_ugL_mean, CO2flux_umolm2s_mean))


# ------ 3. Reforecast Code (The Machine Learning Loop) ------
all_reforecasts <- NULL

message("--- Generating Historical Reforecasts ---")
for(curr_date in reforecast_dates) {
  
  ref_date <- as_date(curr_date)
  curr_site <- focal_sites[1]
  
  message(paste("Processing Reforecast for:", ref_date))
  
  # A. Filter training data to ONLY include days before the reference date
  site_target <- targets_rf |> filter(site_id == curr_site, datetime < ref_date)
  noaa_date <- ref_date - days(1)
  
  # B. Get the historical Stage 2 NOAA forecast for this specific past date
  weather_future_s3 <- tryCatch({
    vera4castHelpers::noaa_stage2(start_date = as.character(noaa_date)) |> 
      dplyr::filter(datetime >= ref_date, site_id == curr_site, variable %in% met_variables) |> 
      dplyr::collect()
  }, error = function(e) { return(data.frame()) })
  
  if(nrow(weather_future_s3) == 0) next 
  
  noaa_future_site <- weather_future_s3 |> 
    mutate(datetime = as_date(datetime)) |> 
    group_by(datetime, site_id, parameter, variable) |> 
    summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
    mutate(prediction = prediction - 273.15) |> pivot_wider(names_from = variable, values_from = prediction) |> 
    select(any_of(c('datetime', 'site_id', met_variables, 'parameter'))) |>
    # Attach covariates to the future drivers
    left_join(model_drivers |> select(datetime, site_id, all_of(covariate_variables)), by = c("datetime", "site_id")) |> 
    tidyr::fill(all_of(covariate_variables), .direction = "down")
  
  if(nrow(site_target) < 10 | nrow(noaa_future_site) < horizon_days) next 
  
  # C. Train the Random Forest on the historical subset
  rf_mod <- parsnip::rand_forest(mode = "regression") |> 
    parsnip::set_engine("ranger") 
  
  rf_recipe <- recipes::recipe(co2 ~ co2_yday + air_temperature + Chla_ugL_mean + CO2flux_umolm2s_mean, data = site_target)
  
  rf_wflow <- workflows::workflow() |> 
    workflows::add_model(rf_mod) |> 
    workflows::add_recipe(rf_recipe)
  
  rf_fit <- rf_wflow |> parsnip::fit(data = site_target)
  
  # D. Setup Ensemble Uncertainty
  forecasted_dates <- seq(from = ref_date, to = max(noaa_future_site$datetime), by = "day")
  curr_co2 <- tail(site_target$co2, 1) 
  
  # Base initial condition uncertainty (Using your 3.7 empirical error)
  ic_uc <- rnorm(n = n_members, mean = curr_co2, sd = 3.7) 
  ic_df <- tibble(forecast_date = rep(ref_date, times = n_members), ensemble_member = 1:n_members, forecast_variable = target_variable, value = ic_uc, uc_type = "total")
  
  forecast_total_unc <- tibble(
    forecast_date = rep(forecasted_dates, times = n_members),
    ensemble_member = rep(1:n_members, each = length(forecasted_dates)),
    forecast_variable = target_variable, value = as.double(NA), uc_type = "total"
  ) |> rows_update(ic_df, by = c("forecast_date","ensemble_member","forecast_variable","uc_type")) 
  
  # E. Step forward through time using the ML model
  for(d in 2:length(forecasted_dates)) {
    temp_pred <- forecast_total_unc |> filter(forecast_date == forecasted_dates[d])
    temp_lag <- forecast_total_unc |> filter(forecast_date == forecasted_dates[d-1])
    temp_driv <- noaa_future_site |> filter(datetime == forecasted_dates[d])
    
    if(nrow(temp_driv) == 0) next
    temp_driv_31 <- temp_driv[rep(1:nrow(temp_driv), length.out = n_members), ]
    
    # Package data for the predict function
    pred_df <- data.frame(
      co2_yday = temp_lag$value,
      air_temperature = temp_driv_31$air_temperature,
      Chla_ugL_mean = temp_driv_31$Chla_ugL_mean,
      CO2flux_umolm2s_mean = temp_driv_31$CO2flux_umolm2s_mean
    )
    
    # Predict and assign
    rf_preds <- predict(rf_fit, new_data = pred_df)
    temp_pred$value <- pmax(rf_preds$.pred, 0)
    
    forecast_total_unc <- forecast_total_unc |> rows_update(temp_pred, by = c("forecast_date","ensemble_member","forecast_variable","uc_type"))
  }
  
  # F. Format for VERA
  curr_reforecast_df <- forecast_total_unc |>
    filter(forecast_date > ref_date) |>
    rename(datetime = forecast_date, parameter = ensemble_member, prediction = value) |>
    mutate(site_id = curr_site, variable = target_variable, reference_datetime = ref_date) |>
    select(datetime, reference_datetime, site_id, parameter, prediction, variable)
  
  all_reforecasts <- dplyr::bind_rows(all_reforecasts, curr_reforecast_df)
}

user_df_EFI <- all_reforecasts |>
  mutate(model_id = my_model_id, family = 'ensemble', duration = 'P1D', project_id = 'vera4cast', parameter = as.character(parameter)) |>
  select(project_id, model_id, datetime, reference_datetime, duration, site_id, family, parameter, variable, prediction)


# ------ 4. Scoring (Model Only) ------
message("--- Scoring Random Forest Model ---")

# Score4cast requires targets to be formatted perfectly to match your forecast
formatted_targets <- targets_co2 |> 
  mutate(variable = target_variable, observation = co2)

# Score just YOUR model (skipping the missing baselines)
scores_df <- score4cast::score(forecast = user_df_EFI, target = formatted_targets)

scoreboard <- scores_df |>
  group_by(model_id) |>
  summarize(mean_crps = round(mean(crps, na.rm = TRUE), 3)) |>
  arrange(mean_crps) 

message("\n=========================================")
message("        2023 REFORECAST SCOREBOARD       ")
message("=========================================")
print(scoreboard |> as.data.frame())
message("=========================================\n")


# ------ 5. Plot CRPS over time ------
message("--- Calculating Mean CRPS by Horizon ---")

horizon_scores <- scores_df |>
  mutate(horizon_days = as.numeric(horizon)) |>
  group_by(model_id, horizon_days) |>
  summarize(mean_crps = mean(crps, na.rm = TRUE), .groups = "drop") |>
  filter(!is.na(horizon_days), horizon_days > 0)

message("--- Generating Plot ---")

horizon_plot <- ggplot(horizon_scores, aes(x = horizon_days, y = mean_crps, color = model_id)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5, alpha = 0.9) +
  theme_bw() +
  scale_color_manual(values = c("dodgerblue4")) + # Just one color for your model
  labs(
    title = "Forecast Performance over Time (FCRE CO2)",
    subtitle = "Mean CRPS by Forecast Horizon (Random Forest)",
    x = "Horizon (Days into the future)",
    y = "Mean CRPS",
    color = "Model"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

print(horizon_plot)