library(data.table)
library(tidyverse)

lifelines_2fl <- fread("data/lifeline/lifeline.2FL.hg38.glm.linear.gz")
# original name format: rsID
# desired name format: CHR_POS_REF_ALT
## effect allele defined as the alternate allele, other allele defined as ref allele
lifelines_2fl$ID <- paste(lifelines_2fl$`#CHROM`, 
                         lifelines_2fl$POS, 
                         lifelines_2fl$REF, 
                         lifelines_2fl$ALT, 
                         sep = "_")

# add sample size column
lifelines_2fl <- lifelines_2fl %>%
  mutate(N = 433)

# check that beta & se columns don't have NAs
summary(lifelines_2fl$BETA)
summary(lifelines_2fl$SE)

### loop for all files
lifeline_files <- list.files(path = "data/lifeline",
                          pattern = "\\hg38.glm.linear.gz$",
                          full.names = TRUE)

clean_lifeline <- function(file, output_file) {
  df <- read_tsv(file)
  
  df$ID <- paste(df$`#CHROM`, 
                 df$POS, 
                 df$REF, 
                 df$ALT, 
                 sep = "_")
  
  df <- df %>% mutate(N = 433)
  
  write_tsv(df, output_file)
  invisible(df)
}

for (i in seq_along(lifeline_files)) {
  clean_lifeline(
    file = lifeline_files[i],
    output_file = sub("\\.glm.linear.gz$", ".cleaned.tsv", lifeline_files[i])
  )
}

# checking output
clean_6SL <- read_tsv("data/lifeline/lifeline.6SL.hg38.cleaned.tsv")
