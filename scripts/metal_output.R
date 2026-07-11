library(data.table)
library(tidyverse)
library(qqman)
library(purrr)

### phenotype venn diagram
library(ggVennDiagram)

data_list <- list(
  MILK = c("2'FL", "3FL", "3'SL", "6'SL", "DFLNH", "DFLNT", "DFLac", "DSLNH", "DSLNT", "FDSLNH", "FLNH", "LNFP I", "LNFP II", "LNFP III", "LNH", "LNT", "LNnT", "LSTb", "LSTc", "Fucosylated Total", "Sialylated Total", "Total HMOs"),
  CHILD = c("2'FL", "3FL", "3'SL", "6'SL", "DFLNH", "DFLNT", "DFLac", "DSLNH", "DSLNT", "FDSLNH", "FLNH", "LNFP I", "LNFP II", "LNFP III", "LNH", "LNT", "LNnT", "LSTb", "LSTc", "Fucosylated Total", "Sialylated Total", "Total HMOs"),
  Lifelines = c("2'FL", "3FL", "3'SL", "6'SL", "DSLNT", "LNFP I", "LNFP II", "LNFP III", "LNH", "LNT", "LNnT", "LSTb", "LSTc", "Fucosylated Total", "Sialylated Total", "Total HMOs")
)

ggVennDiagram(data_list) + 
  scale_fill_gradient(low = "#F4FAFF", high = "#0077B6") +
  labs(title = "Shared HMO Phenotypes Across Cohorts")

### first looking at the SAMPLESIZE output
meta_dat_old <- read.table("results/METAANALYSIS1.TBL", header=TRUE)
meta_dat_old[order(meta_dat_old$P.value), ][1:20, ]

# now looking at STDERR output
meta_dat <- read.table("results/METAANALYSIS1_stderr.TBL", header=TRUE)

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
png(filename = "results/metal_2fl_qqplot.png", width = 800, height = 600, res = 100)
qq(meta_3$P.value)
dev.off()

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

### now looking at all phenotype output
metal_files <- list.files("results/all_cohorts", pattern = "^metal_.*\\.tbl$", full.names = TRUE)
milkchild_files <- list.files("results/milk_child_cohorts", pattern = "^metal_.*\\.tbl$", full.names = TRUE)
milkchild_files <- milkchild_files[-21]

assign_loci <- function(df, window_bp = 500000) {
  df %>%
    arrange(CHR, BP) %>%
    group_by(CHR) %>%
    mutate(
      gap = BP - lag(BP),
      new_locus = if_else(is.na(gap) | gap > window_bp, 1L, 0L),
      locus_id = cumsum(new_locus)
    ) %>%
    ungroup() %>%
    select(-gap, -new_locus)
}

summary_tbl <- function(f, window_bp = 500000) {
  dat <- fread(f)
  
  dat <- dat %>%
    separate(
      MarkerName,
      into = c("CHR", "BP"),
      sep = "_",
      extra = "drop",
      remove = FALSE
    ) %>%
    mutate(
      Phenotype = sub("^metal_|\\.tbl$", "", basename(f)),
      CHR = as.character(CHR),
      BP  = as.numeric(BP),
      P   = as.numeric(`P-value`)
    ) %>%
    filter(!is.na(P), P > 0, P < 5e-8, !is.na(CHR), !is.na(BP))
  
  dat <- assign_loci(dat, window_bp = window_bp)
  
  top_hits <- dat %>%
    group_by(CHR, locus_id) %>%
    slice_min(order_by = P, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      Phenotype = Phenotype,
      CHR = CHR,
      BP = BP,
      SNP = MarkerName,
      Effect_Allele = Allele2,
      Other_Allele = Allele1,
      FreqAvg = Freq1,
      FreqSE = FreqSE,
      Beta = Effect,
      SE = StdErr,
      P = signif(P, 3),
      Direction = Direction,
      HetI2 = HetISq,
      HetPVal = HetPVal
    )
  
  top_hits
}

summ_df <- map_dfr(metal_files, summary_tbl)
summ_df$Phenotype <- sub("1\\.tbl$", "", summ_df$Phenotype)

library(gt)
summ_df %>%
  gt() %>%
  gtsave(filename = "results/lead_snps.html")

### forest plots
forest_plt_df <- function(file_name, obs_id) {
  
  file1 <- paste0("data/milk/milk_", file_name, "_clean.mlma")
  file2 <- paste0("data/child/child_", file_name, "_hg38_cleaned.tsv")
  file3 <- paste0("data/lifeline/lifeline.", file_name, ".hg38.cleaned.tsv")
  
  study1 <- fread(file1) %>%
    filter(SNP == obs_id) %>%
    transmute(
      Phenotype = file_name,
      SNP = SNP,
      Beta = b,
      SE = se,
      Cohort = "MILK"
    )
  
  study2 <- fread(file2) %>%
    filter(rs_id == obs_id) %>%
    transmute(
      Phenotype = file_name,
      SNP = rs_id,
      Beta = beta,
      SE = standard_error,
      Cohort = "CHILD"
    )
  
  study3 <- fread(file3) %>%
    filter(ID == obs_id) %>%
    transmute(
      Phenotype = file_name,
      SNP = ID,
      Beta = BETA,
      SE = SE,
      Cohort = "LIFELINE"
    )
  
  if (nrow(study1) == 0)
    warning(obs_id, " not found in ", file1)
  
  if (nrow(study2) == 0)
    warning(obs_id, " not found in ", file2)
  
  if (nrow(study3) == 0)
    warning(obs_id, " not found in ", file3)
  
  bind_rows(study1, study2, study3)
}

plot_forest_df <- map2_dfr(summ_df$Phenotype, summ_df$Top_SNP, forest_plt_df)
summ_df$Cohort <- "Meta"
plot_forest_df <- rbind(
  plot_forest_df,
  summ_df[, names(plot_forest_df), drop = FALSE]
)

plot_forest_df <- plot_forest_df %>%
  mutate(
    lo = Beta - 1.96 * SE,
    hi = Beta + 1.96 * SE,
    grp = interaction(Phenotype, SNP, drop = TRUE)
  ) %>%
  group_by(grp) %>%
  mutate(Cohort = factor(Cohort, levels = rev(unique(Cohort))),
         Phenotype = droplevels(factor(Phenotype)),
         SNP = droplevels(factor(SNP))) %>%
  ungroup()

ggplot(plot_forest_df, aes(x = Beta, y = Cohort)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
  geom_point(size = 2) +
  facet_wrap(~ Phenotype + SNP, scales = "free_y") +
  labs(x = "Estimate", y = NULL) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold")
  )

### manhattan plot
read_tbl <- function(f) {
  dt <- fread(f)
  
  dt <- dt %>%
    separate(MarkerName,
             into = c("CHR", "BP"),
             sep = "_",
             extra = "drop",
             remove = FALSE) %>%
    mutate(
      trait = sub("^metal_|\\.tbl$", "", basename(f)),
      CHR = as.character(`CHR`),
      BP  = as.numeric(BP),
      P   = as.numeric(`P-value`)
    ) %>%
    filter(!is.na(P), P > 0)
  
  dt
}

df <- map_dfr(metal_files, read_tbl)

df <- df %>%
  mutate(CHR = gsub("^chr", "", CHR)) %>%   # removes "chr" if present
  mutate(CHR = factor(CHR, levels = c(as.character(1:22), "X"))) %>%
  arrange(CHR, BP)

chr_info <- df %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(tot = cumsum(chr_len) - chr_len)

df <- df %>%
  left_join(chr_info, by = "CHR") %>%
  mutate(BP_cum = BP + tot,
         logP = -log10(P))

axis_df <- chr_info %>%
  mutate(center = tot + chr_len / 2)

p <- ggplot(df, aes(x = BP_cum, y = logP, color = trait)) +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
  geom_hline(yintercept = -log10((5e-8)/6), linetype = "dashed") +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR) +
  labs(x = "Chromosome", y = expression(-log[10](P)), color = "Trait") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 3))
ggsave("results/overlay_manhattan_metal_milkchild.png", plot = p,
       width = 14,
       height = 6,
       dpi = 300)

# looking into some loci
sig_df <- df %>%
  dplyr::filter(P <= 5e-8) %>%
  group_by(trait, CHR) %>%
  arrange(P)

het_summ <- sig_df %>%
  summarise(
  n_sig = n(),
  prop_q_p_gt_0.05 = mean(HetPVal > 0.05, na.rm = TRUE),
  mean_HetISq = mean(HetISq, na.rm = TRUE),
  prop_i2_gt_0.5 = mean(HetISq > 50, na.rm = TRUE)
)

summ_df <- df %>%
  dplyr::filter(P <= 5e-8,
         CHR != 19,
         trait != "Sec_stat_child_milk1.tbl") %>%
  distinct(MarkerName, CHR, BP, P, trait, Direction) %>%
  group_by(CHR, trait) %>%
  arrange(P)

library(locuszoomr)
library(ensembldb)
library(EnsDb.Hsapiens.v75)

df_dslnt <- fread("results/all_cohorts/metal_DSLNT1.tbl")
df_dslnt <- df_dslnt %>%
  separate(MarkerName,
           into = c("chrom", "pos"),
           sep = "_",
           extra = "drop",
           remove = FALSE) %>%
  mutate(chrom = as.numeric(chrom),
         pos = as.numeric(pos),
         id = paste0("chr", chrom, ":", pos))
df_dslnt <- df_dslnt %>%
  mutate(`P-value` = as.numeric(`P-value`))
df_dslnt <- df_dslnt[!is.na(df_dslnt$`P-value`), ]

loc <- locuszoomr::locus(data = df_dslnt,
                         seqname = "17",
                         xrange = c(6885193, 7885193),
                         ens_db = "EnsDb.Hsapiens.v75",
                         labs = "id", p = "P-value")
LDlinkR::LDproxy(
  snp = "chr17:7385193",
  pop = "CEU",
  token = "d238732a0b1c",
  genome_build = "grch38"
)
loc$index_snp <- "chr17:7385193"
loc <- locuszoomr::link_LD(
  loc,
  token = "d238732a0b1c",
  method = "proxy",
  genome_build = "grch38"
)

png("results/all_cohorts/locus_plot_dslnt.png", width = 800, height = 600)
locuszoomr::locus_plot(loc)
dev.off()

### chr4 locus
df_lnnt <- fread("results/all_cohorts/metal_LNnT1.tbl")
df_lnnt <- df_lnnt %>%
  separate(MarkerName,
           into = c("chrom", "pos"),
           sep = "_",
           extra = "drop",
           remove = FALSE) %>%
  mutate(chrom = as.numeric(chrom),
         pos = as.numeric(pos),
         id = paste0("chr", chrom, ":", pos))
df_lnnt <- df_lnnt %>%
  mutate(`P-value` = as.numeric(`P-value`))
df_lnnt <- df_lnnt[!is.na(df_lnnt$`P-value`), ]

loc <- locuszoomr::locus(data = df_lnnt,
                         seqname = "4",
                         xrange = c(40826544, 41826544),
                         ens_db = "EnsDb.Hsapiens.v75",
                         labs = "id", p = "P-value")
LDlinkR::LDproxy(
  snp = "chr4:41326544",
  pop = "CEU",
  token = "d238732a0b1c",
  genome_build = "grch38"
)
loc$index_snp <- "chr4:41326544"
loc <- locuszoomr::link_LD(
  loc,
  token = "d238732a0b1c",
  method = "proxy",
  genome_build = "grch38"
)

png("results/all_cohorts/locus_plot_lnnt.png", width = 800, height = 600)
locuszoomr::locus_plot(loc)
dev.off()

### chr6 locus
df_dslnt <- fread("results/all_cohorts/metal_DSLNT1.tbl")
df_dslnt <- df_dslnt %>%
  separate(MarkerName,
           into = c("chrom", "pos"),
           sep = "_",
           extra = "drop",
           remove = FALSE) %>%
  mutate(chrom = as.numeric(chrom),
         pos = as.numeric(pos),
         id = paste0("chr", chrom, ":", pos))
df_dslnt <- df_dslnt %>%
  mutate(`P-value` = as.numeric(`P-value`))
df_dslnt <- df_dslnt[!is.na(df_dslnt$`P-value`), ]

loc <- locuszoomr::locus(data = df_dslnt,
                         seqname = "6",
                         xrange = c(112724735, 113724735),
                         ens_db = "EnsDb.Hsapiens.v75",
                         labs = "id", p = "P-value")
LDlinkR::LDproxy(
  snp = "chr6:113224735",
  pop = "CEU",
  token = "d238732a0b1c",
  genome_build = "grch38"
)
loc$index_snp <- "chr17:7385193"
loc <- locuszoomr::link_LD(
  loc,
  token = "d238732a0b1c",
  method = "proxy",
  genome_build = "grch38"
)

png("results/all_cohorts/locus_plot_dslnt_chr6.png", width = 800, height = 600)
locuszoomr::locus_plot(loc)
dev.off()
