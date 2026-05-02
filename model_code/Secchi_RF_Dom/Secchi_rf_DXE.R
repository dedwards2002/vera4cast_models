# 1. Load the required packages
library(tidyverse)
library(lubridate)
library(vera4castHelpers)

# 2. Set arguments
my_forecast_date <- Sys.Date() - lubridate::days(1)
my_model_id <- 'Dom_RF_Secchi' 

# 3. Source the math function from your R folder
source('./R/Secchi_Dom/Secchi_RF_Model_Code.R')

# 4. Run the function!
final_forecast <- Dom_RF_Secchi_model_function(forecast_date = my_forecast_date, model_id = my_model_id)

# ---------------------------------------------------------
# SAVE AND SUBMIT BLOCK
# ---------------------------------------------------------

# Step A: Guarantee the folder exists BEFORE trying to save anything
dir.create("model_output/secchi_DXE", showWarnings = FALSE, recursive = TRUE)

# Step B: Save the file using the clean path (no leading ./)
write.csv(final_forecast, "model_output/secchi_DXE/Secchi_rf_df.csv", row.names = FALSE)

# Step C: Validate the file using the clean path
vera4castHelpers::forecast_output_validator("model_output/secchi_DXE/Secchi_rf_df.csv")

# Step D: Submit the file using the clean path
vera4castHelpers::submit("model_output/secchi_DXE/Secchi_rf_df.csv", s3_region = "submit", s3_endpoint = "ltreb-reservoirs.org", first_submission = FALSE)