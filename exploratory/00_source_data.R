# Load required libraries
library(googledrive)
library(fs)

# Set up the Google Drive link
drive_download("https://drive.google.com/uc?id=1uRoPmYxQG2M7F5Q8cDfCmG6yXYnGkR4d/", 
               path = "data/raw/dataset.csv",
               overwrite = TRUE)
