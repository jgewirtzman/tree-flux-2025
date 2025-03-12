# Load required libraries
library(googledrive)
library(fs)

# Set up the Google Drive link
drive_download("https://drive.google.com/uc?id=1vowUAjBXJ26Y3nlrYsL-2bqjILww1lrF/", 
               path = "data/raw/flux_dataset.csv",
               overwrite = TRUE)


drive_download("https://drive.google.com/uc?id=1vW8xBF5BuHsr63mRxoo1CZntb388llNi/", 
               path = "data/raw/ymf_dataset.csv",
               overwrite = TRUE)

drive_download("https://drive.google.com/uc?id=1wYVEZEJaI1PsTbuo7SZDFzH2Hypolzdc/", 
               path = "data/raw/ymf_dataset.csv",
               overwrite = TRUE)

drive_download("https://drive.google.com/uc?id=1HvdNJqNjvzyRYNO4Wv1_kla8TlMUsTGI/", 
               path = "data/raw/wtd_ems-met.csv",
               overwrite = TRUE)

