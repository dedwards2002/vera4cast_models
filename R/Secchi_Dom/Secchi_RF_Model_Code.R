Dom_RF_Secchi_model_function <- function(forecast_date, model_id) {
  
  forecast_start_date <- forecast_date
  
  # --- 1. Variables & Setup ---
  focal_sites <- c("fcre")
  target_variable <- "Secchi_m_sample"
  
  covariate_variables <- c("SchmidtStability_Jm2_mean", "Temp_C_mean") 
  met_variables <- c("air_temperature", "precipitation_flux", "eastward_wind", "northward_wind")
  n_members <- 31
  forecast_df <- NULL
  
  message("--- Downloading Target Data ---")
  targets_url <- "https://amnh1.osn.mghpcc.org/bio230121-bucket01/vera4cast/targets/project_id=vera4cast/duration=P1D/daily-insitu-targets.csv.gz"
  targets <- read_csv(targets_url, show_col_types = FALSE) |> 
    filter(site_id %in% focal_sites) |> mutate(datetime = as_date(datetime)) 
  
  message("--- Downloading NOAA Weather Data ---")
  weather_past <- vera4castHelpers::noaa_stage3() |> 
    dplyr::filter(site_id %in% focal_sites, datetime >= ymd('2018-01-01'), variable %in% met_variables) |> 
    dplyr::collect()
  
  weather_past_daily <- weather_past |> 
    mutate(datetime = as_date(datetime)) |> group_by(datetime, site_id, variable) |> 
    summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
    mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
    pivot_wider(names_from = variable, values_from = prediction)
  
  weather_future <- vera4castHelpers::noaa_stage2(start_date = as.character(forecast_date - days(1))) |> 
    dplyr::filter(datetime >= forecast_date, site_id %in% focal_sites, variable %in% met_variables) |> collect()
  
  if(nrow(weather_future) == 0) {
    weather_future <- vera4castHelpers::noaa_stage2(start_date = as.character(forecast_date - days(2))) |> 
      dplyr::filter(datetime >= forecast_date, site_id %in% focal_sites, variable %in% met_variables) |> collect()
  }
  
  weather_future_daily <- weather_future |> mutate(datetime = as_date(datetime)) |> 
    group_by(datetime, site_id, parameter, variable) |> summarize(prediction = mean(prediction, na.rm = TRUE), .groups = "drop") |> 
    mutate(prediction = ifelse(variable == "air_temperature", prediction - 273.15, prediction)) |> 
    pivot_wider(names_from = variable, values_from = prediction) |> 
    select(any_of(c('datetime', 'site_id', met_variables, 'parameter')))
  
  message("--- Cleaning and Merging ---")
  # Use across to summarize all met variables properly 
  weather_combined <- bind_rows(weather_past_daily, weather_future_daily |> select(datetime, site_id, any_of(met_variables)) |> group_by(datetime, site_id) |> summarize(across(any_of(met_variables), ~mean(.x, na.rm = TRUE)), .groups = 'drop')) |> distinct(datetime, site_id, .keep_all = TRUE) |> arrange(site_id, datetime)
  
  covariates_historical <- targets |> filter(variable %in% covariate_variables) |> select(datetime, site_id, variable, observation) |> group_by(datetime, site_id, variable) |> summarize(observation = mean(observation, na.rm = TRUE), .groups = "drop") |> pivot_wider(names_from = variable, values_from = observation)
  
  model_drivers <- weather_combined |> left_join(covariates_historical, by = c("datetime", "site_id")) |> distinct(datetime, site_id, .keep_all = TRUE) |> group_by(site_id) |> arrange(datetime) |> tidyr::fill(all_of(covariate_variables), .direction = "down") |> ungroup()
  
  # Removed depth filter for Secchi
  targets_secchi <- targets |> filter(variable == target_variable) |> group_by(datetime, site_id) |> summarize(secchi = mean(observation, na.rm = TRUE), .groups = "drop")
  
  # Ensure complete.cases checks new met variables
  targets_lm <- targets_secchi |> left_join(model_drivers, by = c("datetime", "site_id")) |> arrange(site_id, datetime) |> group_by(site_id) |> mutate(secchi_yday = lag(secchi, 1)) |> ungroup() |> filter(complete.cases(secchi, secchi_yday, air_temperature, precipitation_flux, eastward_wind, northward_wind, SchmidtStability_Jm2_mean, Temp_C_mean))
  
  # --- 2. Fit model & generate forecast ---
  message("--- Generating Forecast Ensembles ---")
  for(i in 1:length(focal_sites)) {  
    
    curr_site <- focal_sites[i]
    site_target <- targets_lm |> filter(site_id == curr_site)
    
    noaa_future_site <- weather_future_daily |> filter(site_id == curr_site) |> left_join(model_drivers |> select(datetime, site_id, all_of(covariate_variables)), by = c("datetime", "site_id")) |> distinct(datetime, site_id, parameter, .keep_all = TRUE)
    
    # Train the Random Forest
    rf_mod <- parsnip::rand_forest(mode = "regression") |> 
      parsnip::set_engine("ranger") 
    
    rf_recipe <- recipes::recipe(secchi ~ secchi_yday + air_temperature + precipitation_flux + eastward_wind + northward_wind + SchmidtStability_Jm2_mean + Temp_C_mean, data = site_target)
    
    rf_wflow <- workflows::workflow() |> 
      workflows::add_model(rf_mod) |> 
      workflows::add_recipe(rf_recipe)
    
    rf_fit <- rf_wflow |> 
      parsnip::fit(data = site_target)
    
    forecasted_dates <- seq(from = forecast_date, to = max(noaa_future_site$datetime), by = "day")
    curr_secchi <- tail(site_target$secchi, 1) 
    
    # Initial Condition Uncertainty (Adjusted for Secchi depth scale)
    ic_uc <- rnorm(n = n_members, mean = curr_secchi, sd = 0.5) 
    ic_uc <- pmax(ic_uc, 0.1) # Prevent negative starting secchi depths
    
    ic_df <- tibble(forecast_date = rep(forecast_date, times = n_members), ensemble_member = 1:n_members, forecast_variable = target_variable, value = ic_uc, uc_type = "total")
    
    forecast_total_unc <- tibble(forecast_date = rep(forecasted_dates, times = n_members), ensemble_member = rep(1:n_members, each = length(forecasted_dates)), forecast_variable = target_variable, value = as.double(NA), uc_type = "total") |> rows_update(ic_df, by = c("forecast_date","ensemble_member","forecast_variable","uc_type")) 
    
    # The Machine Learning Forecast Loop
    for(d in 2:length(forecasted_dates)) {
      temp_pred <- forecast_total_unc |> filter(forecast_date == forecasted_dates[d])
      temp_lag <- forecast_total_unc |> filter(forecast_date == forecasted_dates[d-1])
      temp_driv <- noaa_future_site |> filter(datetime == forecasted_dates[d])
      
      temp_driv_31 <- temp_driv[rep(1:nrow(temp_driv), length.out = n_members), ]
      
      pred_df <- data.frame(
        secchi_yday = temp_lag$value,
        air_temperature = temp_driv_31$air_temperature,
        precipitation_flux = temp_driv_31$precipitation_flux,
        eastward_wind = temp_driv_31$eastward_wind,
        northward_wind = temp_driv_31$northward_wind,
        SchmidtStability_Jm2_mean = temp_driv_31$SchmidtStability_Jm2_mean,
        Temp_C_mean = temp_driv_31$Temp_C_mean
      )
      
      # Generate predictions
      rf_preds <- predict(rf_fit, new_data = pred_df)
      
      # Assign predictions, preventing negative Secchi depth
      temp_pred$value <- pmax(rf_preds$.pred, 0.1)
      
      forecast_total_unc <- forecast_total_unc |> rows_update(temp_pred, by = c("forecast_date","ensemble_member","forecast_variable","uc_type"))
    }
    
    # depth_m to be NA
    curr_site_df <- forecast_total_unc |> filter(forecast_date > forecast_start_date) |> rename(datetime = forecast_date, parameter = ensemble_member, prediction = value) |> mutate(site_id = curr_site, variable = target_variable, depth_m = as.numeric(NA)) |> select(datetime, site_id, depth_m, parameter, prediction, variable)
    forecast_df <- dplyr::bind_rows(forecast_df, curr_site_df)
  }
  
  # --- 3. Format to VERA standard ---
  forecast_df_EFI <- forecast_df %>%
    mutate(model_id = model_id, reference_datetime = forecast_date, family = 'ensemble', duration = 'P1D', parameter = as.character(parameter), project_id = 'vera4cast') %>%
    select(project_id, model_id, datetime, reference_datetime, duration, site_id, depth_m, family, parameter, variable, prediction)
  
  return(forecast_df_EFI)
}