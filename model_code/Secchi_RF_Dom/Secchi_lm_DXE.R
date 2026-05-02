# 1. Load the required packages
library(tidyverse)
library(lubridate)
library(vera4castHelpers)

# 2. Set arguments
my_forecast_date <- Sys.Date() - lubridate::days(1)
my_model_id <- 'Secchi_lm_Dom'

# 3. Source the math function from  R folder
source('./R/Secchi_Dom/Secchi_LM_model_code.R')

# 4. Run the function!
final_forecast <- Dom_LM_Secchi_model_function(forecast_date = my_forecast_date, model_id = my_model_id)

# ---------------------------------------------------------
# SAVE AND SUBMIT BLOCK
# ---------------------------------------------------------

# Step A: Guarantee the folder exists before trying to save anything
dir.create("model_output/secchi_DXE", showWarnings = FALSE, recursive = TRUE)

# Step B: Save the file
write.csv(final_forecast, "model_output/secchi_DXE/Secchi_df.csv", row.names = FALSE)

# Step C: Validate the file 
vera4castHelpers::forecast_output_validator("model_output/secchi_DXE/Secchi_df.csv")

# Step D: Submit the file 
vera4castHelpers::submit("model_output/secchi_DXE/Secchi_df.csv", s3_region = "submit", s3_endpoint = "ltreb-reservoirs.org", first_submission = FALSE)