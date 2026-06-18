library(data.table)
library(tidyverse)
library(qqman)

meta_dat <- read.table("results/METAANALYSIS1.TBL", header=TRUE)
meta_dat[order(meta_dat$P.value), ][1:20, ]

# looks like MILK study contributely solely to many top hits
# filter to hits where all 3 studies contributed
meta_3 <- meta_dat %>%
  mutate(P.value = as.numeric(P.value)) %>%
  filter(!str_detect(Direction, fixed("?"))) %>%
  arrange(`P.value`)

# how much directional concordance is there
meta_3 %>%
  filter(`P.value` < 1e-3) %>%
  summarise(
    n_total = n(),
    n_concordant = sum(Direction %in% c("+++", "---")),
    concordance_prop = n_concordant / n_total
  )

# qq plot
plot <- qq(meta_3$P.value)

# manhattan plot
meta_3 <- meta_3 %>%
  separate(MarkerName,
           into = c("CHR", "BP"),
           sep = "_",
           extra = "drop",
           remove = FALSE) %>%
  mutate(CHR = as.numeric(CHR),
         BP = as.numeric(BP))

png("results/metal1_manhattan.png", width = 480*4, height = 480*2, pointsize = 24)
manhattan(meta_3, 
          chr = "CHR",
          bp = "BP",
          p = "P.value",
          snp = "MarkerName", 
          genomewideline = -log10(5e-8),
          suggestiveline = -log10(1e-5),
          annotatePval = 5*10**-8 )
dev.off()
