library(tidyverse)
library(data.table)

# first examine one file
child_2fl <- read_tsv("data/child/child_2FL_hg38.tsv")
#remove unmapped variants
unmapped_2fl <- read_tsv(
  "data/child/2FL.unmapped.bed",
  col_names = FALSE,
  comment = "#" 
) 
colnames(unmapped_2fl) <- c("chr", "start", "end", "rs_id")

child_2fl <- child_2fl %>%
  filter(!(rs_id %in% unmapped_2fl$rs_id))

# original name format: CHR:POS:ALT:REF using hg19 coordinates
# desired name format: CHR_POS_REF_ALT using hg38 coordinates
# effect allele defined as the alternate allele, other allele defined as ref allele
# remove "chr" from chr names
child_2fl <- child_2fl %>%
  filter(!grepl("_alt|_random", chromosome)) %>%
  mutate(chromosome = sub("^chr", "", chromosome))
child_2fl$rs_id <- paste(child_2fl$chromosome, 
                         child_2fl$base_pair_location, 
                         child_2fl$other_allele, 
                         child_2fl$effect_allele, 
                         sep = "_")

# check that beta & se columns don't have NAs
summary(child_2fl$beta)
summary(child_2fl$standard_error)

### loop for all files
child_files <- list.files(path = "data/child",
                          pattern = "\\hg38.tsv$",
                          full.names = TRUE)

child_unmapped <- list.files(path = "data/child",
                             pattern = "\\.unmapped.bed$",
                             full.names = TRUE)

clean_child <- function(file, unmapped, output_file) {
  df <- read_tsv(file)
  unmapped_df <- read_tsv(
    unmapped,
    col_names = FALSE,
    comment = "#" 
  ) 
  colnames(unmapped_df) <- c("chr", "start", "end", "rs_id")
  
  df <- df %>%
    filter(!(rs_id %in% unmapped_df$rs_id))
  
  df <- df %>%
    filter(!grepl("_alt|_random", chromosome)) %>%
    mutate(chromosome = sub("^chr", "", chromosome))
  
  df$rs_id<- paste(df$chromosome, 
                   df$base_pair_location, 
                   df$other_allele, 
                   df$effect_allele, 
                   sep = "_")
  
  write_tsv(df, output_file)
  invisible(df)
}

for (i in seq_along(child_files)) {
  clean_child(
    file = child_files[i],
    unmapped = child_unmapped[i],
    output_file = sub("\\.tsv$", "_cleaned.tsv", child_files[i])
  )
}

# checking output
clean_LNT <- read_tsv("data/child/child_LNT_hg38_cleaned.tsv")
