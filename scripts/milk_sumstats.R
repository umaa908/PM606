library(tidyverse)
library(data.table)

# first examine one file
milk_2fl <- read_tsv("data/milk/gcta_int_2FL.mlma")

# check that beta & se & p-val columns don't have NAs
summary(milk_2fl$b)
summary(milk_2fl$se)
summary(milk_2fl$p)

milk_files <- list.files(path = "data/milk",
                          pattern = ".mlma$",
                          full.names = TRUE)

clean_milk <- function(file, output_file) {
  df <- read_tsv(file)
  
  df <- df %>%
    filter(!is.na(b),
           !is.na(p),
           str_count(SNP, "_") == 3)
  
  df <- df %>% mutate(N = 349)
  
  write_tsv(df, output_file)
  invisible(df)
}

for (i in seq_along(milk_files)) {
  clean_milk(
    file = milk_files[i],
    output_file = sub("^gcta_int_(.*)\\.mlma$", "milk_\\1_clean.mlma", milk_files[i])
  )
}
