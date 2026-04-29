# ==============================================================================
# SECCHI REFORECAST: Random Forest vs Linear Model vs VERA Baselines (2024)
# ==============================================================================

# ------ Load packages -----
library(tidyverse)
library(lubridate)
library(vera4castHelpers)
library(parsnip)
library(recipes)
library(workflows)
library(score4cast)
library(arrow)

# --- Configuration ---
rf_model_id <- 'Dom_RF_Secchi'
lm_model_id <- 'Dom_LM_Secchi'
focal_sites <- c("fcre")
target_variable <- "Secchi_m_sample"

# Secchi relies on specific stability and temperature covariates
covariate_variables <- c("SchmidtStability_Jm2_mean", "Temp_C_mean") 
met_variables <- c("air_temperature", "precipitation_flux", "eastward_wind", "northward_wind")

reforecast_dates <- seq(ymd('2024-04-01'), ymd('2024-10-31'), by = "week")
horizon_days <- 30 
n_members <- 31

#------- 1. Read Target Data --------
message("--- Downloading VERA Targets ---")
targets_url <- "https://amnh1.osn.mghpcc.org/bio230121-bucket01/vera4cast/targets/project_id=vera4cast/duration=P1D/daily-insitu-targets.csv.gz"
targets <- read_csv(targets_url, show_col_types = FALSE) |>
  filter(site_id %in% focal_sites) |> 
  mutate(datetime = as_date(datetime))

# ------ 2. Download Historical Weather & Drivers ------
message("--- Downloading Historical NOAA Stage 3 ---")
weather_past_s3 <- vera4castHelpers::noaa_stage3() |> 
  dplyr::filter(site_id %in% focal_sites, datetime >= ymd('2018-01-01'), variable %in% met_variables) |> 
  dplyr::collect()

weather_past_daily <- weather_past_s3 |> 
  mutate(datetime = as_date(datetime)) |> 
  group_by(datetime, site_id, variable) |> 
  summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
  mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
  pivot_wider(names_from = variable, values_from = prediction) |> 
  arrange(site_id, datetime)

message("--- Cleaning and Merging Drivers ---")
covariates_historical <- targets |> 
  filter(variable %in% covariate_variables) |> 
  select(datetime, site_id, variable, observation) |> 
  group_by(datetime, site_id, variable) |> 
  summarize(observation = mean(observation, na.rm = TRUE), .groups = "drop") |> 
  pivot_wider(names_from = variable, values_from = observation)

model_drivers <- weather_past_daily |> 
  left_join(covariates_historical, by = c("datetime", "site_id")) |> 
  distinct(datetime, site_id, .keep_all = TRUE) |> 
  group_by(site_id) |> arrange(datetime) |> 
  tidyr::fill(all_of(c(covariate_variables, met_variables)), .direction = "down") |> 
  ungroup()

targets_secchi <- targets |> 
  filter(variable == target_variable) |> 
  group_by(datetime, site_id) |> 
  summarize(secchi = mean(observation, na.rm = TRUE), .groups = "drop")

# Final merged dataset for both LM and RF
targets_lm <- targets_secchi |> 
  left_join(model_drivers, by = c("datetime", "site_id")) |> 
  arrange(site_id, datetime) |> 
  group_by(site_id) |> 
  mutate(secchi_yday = lag(secchi, 1)) |> 
  ungroup() |> 
  filter(complete.cases(secchi, secchi_yday, air_temperature, precipitation_flux, eastward_wind, northward_wind, SchmidtStability_Jm2_mean, Temp_C_mean))

# ==============================================================================
# 3. Reforecast Code (Dual Models)
# ==============================================================================
all_reforecasts <- NULL
curr_site <- focal_sites[1]

for(curr_date in reforecast_dates) {
  
  ref_date <- as_date(curr_date)
  message(paste("Running LM and RF ensembles for reference date:", ref_date))
  
  site_target_train <- targets_lm |> filter(site_id == curr_site, datetime < ref_date)
  if(nrow(site_target_train) < 30) next 
  
  # Fetch future NOAA weather
  forecasted_dates <- seq(from = ref_date, to = ref_date + days(horizon_days), by = "day")
  
  noaa_future_s2 <- tryCatch({
    vera4castHelpers::noaa_stage2(start_date = as.character(ref_date - days(1))) |> 
      dplyr::filter(datetime >= ref_date, site_id == curr_site, variable %in% met_variables) |> collect()
  }, error = function(e) { return(data.frame()) })
  
  if(nrow(noaa_future_s2) == 0) next
  
  noaa_future_site <- noaa_future_s2 |> 
    mutate(datetime = as_date(datetime)) |> 
    group_by(datetime, site_id, parameter, variable) |> 
    summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
    mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
    pivot_wider(names_from = variable, values_from = prediction) |> 
    left_join(model_drivers |> select(datetime, site_id, all_of(covariate_variables)), by = c("datetime", "site_id")) |> 
    distinct(datetime, site_id, parameter, .keep_all = TRUE)
  
  # ------ Train Both Models ------
  # Linear Model
  lm_fit <- lm(secchi ~ secchi_yday + air_temperature + precipitation_flux + eastward_wind + northward_wind + SchmidtStability_Jm2_mean + Temp_C_mean, data = site_target_train)
  lm_sigma <- sd(lm_fit$residuals, na.rm = TRUE)
  
  # Random Forest
  rf_mod <- parsnip::rand_forest(mode = "regression") |> parsnip::set_engine("ranger") 
  rf_recipe <- recipes::recipe(secchi ~ secchi_yday + air_temperature + precipitation_flux + eastward_wind + northward_wind + SchmidtStability_Jm2_mean + Temp_C_mean, data = site_target_train)
  rf_wflow <- workflows::workflow() |> workflows::add_model(rf_mod) |> workflows::add_recipe(rf_recipe)
  rf_fit <- rf_wflow |> parsnip::fit(data = site_target_train)
  
  # ------ Setup Initial Conditions ------
  curr_secchi <- tail(site_target_train$secchi, 1) 
  ic_uc <- rnorm(n = n_members, mean = curr_secchi, sd = 0.5) 
  ic_uc <- pmax(ic_uc, 0.1) 
  
  ic_df <- tibble(forecast_date = rep(ref_date, times = n_members), ensemble_member = 1:n_members, forecast_variable = target_variable, value = ic_uc, uc_type = "total")
  
  lm_forecast_unc <- tibble(forecast_date = rep(forecasted_dates, times = n_members), ensemble_member = rep(1:n_members, each = length(forecasted_dates)), forecast_variable = target_variable, value = as.double(NA), uc_type = "total") |> rows_update(ic_df, by = c("forecast_date","ensemble_member","forecast_variable","uc_type")) 
  
  rf_forecast_unc <- lm_forecast_unc # Duplicate for RF
  
  # ------ Daily Forecast Loop ------
  for(d in 2:length(forecasted_dates)) {
    lm_temp_pred <- lm_forecast_unc |> filter(forecast_date == forecasted_dates[d])
    rf_temp_pred <- rf_forecast_unc |> filter(forecast_date == forecasted_dates[d])
    
    lm_temp_lag <- lm_forecast_unc |> filter(forecast_date == forecasted_dates[d-1])
    rf_temp_lag <- rf_forecast_unc |> filter(forecast_date == forecasted_dates[d-1])
    
    temp_driv <- noaa_future_site |> filter(datetime == forecasted_dates[d])
    if(nrow(temp_driv) == 0) temp_driv <- noaa_future_site |> filter(datetime == forecasted_dates[d-1])
    temp_driv_expanded <- temp_driv[rep(1:nrow(temp_driv), length.out = n_members), ]
    
    # LM Predict
    lm_pred_df <- data.frame(secchi_yday = lm_temp_lag$value, air_temperature = temp_driv_expanded$air_temperature, precipitation_flux = temp_driv_expanded$precipitation_flux, eastward_wind = temp_driv_expanded$eastward_wind, northward_wind = temp_driv_expanded$northward_wind, SchmidtStability_Jm2_mean = temp_driv_expanded$SchmidtStability_Jm2_mean, Temp_C_mean = temp_driv_expanded$Temp_C_mean)
    lm_preds <- predict(lm_fit, newdata = lm_pred_df) + rnorm(n_members, 0, lm_sigma)
    lm_temp_pred$value <- pmax(lm_preds, 0.1)
    
    # RF Predict
    rf_pred_df <- data.frame(secchi_yday = rf_temp_lag$value, air_temperature = temp_driv_expanded$air_temperature, precipitation_flux = temp_driv_expanded$precipitation_flux, eastward_wind = temp_driv_expanded$eastward_wind, northward_wind = temp_driv_expanded$northward_wind, SchmidtStability_Jm2_mean = temp_driv_expanded$SchmidtStability_Jm2_mean, Temp_C_mean = temp_driv_expanded$Temp_C_mean)
    rf_preds <- predict(rf_fit, new_data = rf_pred_df)
    rf_temp_pred$value <- pmax(rf_preds$.pred, 0.1)
    
    # Update DataFrames
    lm_forecast_unc <- lm_forecast_unc |> rows_update(lm_temp_pred, by = c("forecast_date","ensemble_member","forecast_variable","uc_type"))
    rf_forecast_unc <- rf_forecast_unc |> rows_update(rf_temp_pred, by = c("forecast_date","ensemble_member","forecast_variable","uc_type"))
  }
  
  # Clean and tag models
  lm_clean <- lm_forecast_unc |> filter(forecast_date > ref_date) |> mutate(model_id = lm_model_id)
  rf_clean <- rf_forecast_unc |> filter(forecast_date > ref_date) |> mutate(model_id = rf_model_id)
  
  curr_reforecast_df <- bind_rows(lm_clean, rf_clean) |>
    rename(datetime = forecast_date, parameter = ensemble_member, prediction = value) |>
    mutate(site_id = curr_site, variable = target_variable, reference_datetime = ref_date, depth_m = NA) |>
    select(model_id, datetime, reference_datetime, site_id, depth_m, parameter, prediction, variable)
  
  all_reforecasts <- dplyr::bind_rows(all_reforecasts, curr_reforecast_df)
}

# --- Standard Format ---
user_df_EFI <- all_reforecasts |>
  mutate(family = 'ensemble', duration = 'P1D', project_id = 'vera4cast', parameter = as.character(parameter)) |>
  select(project_id, model_id, datetime, reference_datetime, duration, site_id, depth_m, family, parameter, variable, prediction)

# ==============================================================================
# 4. Scoring & Official Baselines
# ==============================================================================

# 1. Score Local Models
targets_score <- targets |> filter(variable == target_variable, site_id == "fcre")
local_scores <- score4cast::score(forecast = user_df_EFI, target = targets_score)

# Prep Local Scores 
local_prep <- local_scores |>
  mutate(horizon_days = as.numeric(horizon), model_id = as.character(model_id)) |>
  select(model_id, horizon_days, crps)

# 2. Fetch VERA Official Baselines via Arrow
message("\n--- Fetching Official VERA Baselines ---")
s3_scores <- arrow::s3_bucket(
  "bio230121-bucket01/vera4cast/scores/parquet", 
  endpoint_override = "amnh1.osn.mghpcc.org", 
  anonymous = TRUE
)

baseline_raw <- arrow::open_dataset(s3_scores) |>
  dplyr::filter(
    variable == target_variable,
    site_id == "fcre",
    model_id %in% c("climatology", "persistenceRW")
  ) |>
  dplyr::collect() |>
  mutate(ref_date_clean = as_date(reference_datetime)) |>
  filter(ref_date_clean %in% as_date(reforecast_dates))

# Prep Baselines by calculating horizon manually 
baseline_prep <- baseline_raw |>
  mutate(
    horizon_days = as.numeric(as_date(datetime) - as_date(reference_datetime)),
    model_id = as.character(model_id)
  ) |>
  select(model_id, horizon_days, crps)

# 3. Combine and Aggregate
horizon_plot_data <- bind_rows(local_prep, baseline_prep) |>
  filter(!is.na(horizon_days), horizon_days > 0) |>
  group_by(model_id, horizon_days) |>
  summarize(mean_crps = mean(crps, na.rm = TRUE), .groups = "drop") |>
  filter(!is.nan(mean_crps))

# Scoreboard
scoreboard <- horizon_plot_data |> 
  group_by(model_id) |> 
  summarize(overall_mean_crps = round(mean(mean_crps, na.rm=TRUE), 3)) |> 
  arrange(overall_mean_crps)

message("\n=========================================")
message("       2024 REFORECAST SCOREBOARD      ")
message("=========================================")
print(scoreboard |> as.data.frame())
message("=========================================\n")

# ==============================================================================
# 5. Plot CRPS over time by model
# ==============================================================================
message("--- Generating Four-Model Comparison Plot ---")

horizon_plot <- ggplot(horizon_plot_data, aes(x = horizon_days, y = mean_crps, color = model_id, group = model_id)) +
  # Points indicate actual sampling days (Secchi is sparse)
  geom_point(size = 2, alpha = 0.5) +
  # LOESS trendline interpolates between sampling dates 
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.5, span = 0.8) +
  theme_bw() +
  scale_color_manual(values = c(
    "climatology" = "gray50", 
    "persistenceRW" = "darkorange", 
    "Dom_LM_Secchi" = "firebrick",
    "Dom_RF_Secchi" = "dodgerblue4"
  )) + 
  labs(
    title = "Forecast Performance over Time (FCRE Secchi Depth)",
    subtitle = "2024 Reforecast: LM vs RF vs VERA Baselines",
    x = "Horizon (Days into the future)",
    y = "Mean CRPS (Meters - Lower is Better)",
    color = "Model"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

print(horizon_plot)