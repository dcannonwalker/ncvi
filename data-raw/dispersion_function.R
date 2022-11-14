## code to prepare `dispersion_function` dataset goes here
load('../RiboSeq_share 2/fits/edger_dpi8_n04.rds')
dispersion_function <- data.frame(dispersion = edger_dpi8_n04$y$tagwise.dispersion,
                 AveLogCPM = edger_dpi8_n04$y$AveLogCPM)
usethis::use_data(dispersion_function, overwrite = TRUE)
