
# 1. Load the required packages
library(tidyverse)
library(lubridate)
library(vera4castHelpers)

# 2. Set arguments
my_forecast_date <- Sys.Date() - lubridate::days(1)
my_model_id <- 'Secchi_rf_Dom'

# 3. Source the math function from your R folder
source('./R/Secchi_Dom/Secchi_RF_Model_Code.R')

# 4. Run the function!
final_forecast <- Dom_RF_Secchi_model_function(forecast_date = my_forecast_date, model_id = my_model_id)

# save and submit
write.csv(final_forecast, './model_output/secchi_DXE/Secchi_rf_df.csv')

vera4castHelpers::forecast_output_validator('./model_output/secchi_DXE/Secchi_rf_df.csv')

vera4castHelpers::submit('./model_output/secchi_DXE/Secchi_rf_df.csv', s3_region = "submit", s3_endpoint = "ltreb-reservoirs.org", first_submission = FALSE)