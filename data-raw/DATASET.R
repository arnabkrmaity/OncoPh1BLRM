## code to prepare `western.trial` dataset goes here

western.trial <- list(DosesAdm = c(1, 2, 4, 8, 15, 30),
                          Npat     = c(5, 6, 5, 9,  8,  4),
                           Ntox     = c(0, 0, 0, 0,  1,  3),
                          Study      = c(1, 1, 1, 1,  1,  1),
                       StudyLabel = c("Western"))

usethis::use_data(western.trial, overwrite = TRUE)



