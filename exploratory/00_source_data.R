# Load required libraries
library(googledrive)
library(fs)

# Set up the Google Drive link
drive_download("https://drive.google.com/uc?id=1vowUAjBXJ26Y3nlrYsL-2bqjILww1lrF/", 
               path = "data/raw/flux_dataset.csv",
               overwrite = TRUE)


