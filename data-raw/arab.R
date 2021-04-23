## code to prepare `arab` dataset goes here

## would like to use function from threadr package,
## but it had an exception so I just copied the function
read_rds_remote <- function(file_remote, quiet = TRUE) {

  # Remote read
  if (grepl("^http", file_remote, ignore.case = TRUE)) {

    # Make temp local copy
    file_local <- basename(file_remote)
    file_local <- file.path(tempdir(), file_local)

    # Download file
    download.file(file_remote, file_local, quiet = quiet, mode = "wb")

    # Load file
    df <- readRDS(file_local)

  } else {

    # Direct read
    df <- readRDS(file_remote)

  }

  # Return
  df

}

arab <- read_rds_remote("http://bioinf.wehi.edu.au/edgeR/UserGuideData/arab.rds")

usethis::use_data(arab, overwrite = TRUE)
