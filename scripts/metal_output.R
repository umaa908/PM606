library(data.table)
library(tidyverse)
library(qqman)
library(purrr)

### phenotype euler diagram
library(ggforce)

# Outer set: MILK = CHILD, 22 elements
# Inner set: LIFELINE, 16 of those 22 elements
# Therefore 6 elements are in A/B but not C

circles <- data.frame(
  x0 = c(0, 0),
  y0 = c(0, 0),
  r  = c(2.2, 1.35),
  set = c("A = B (22)", "C (16)")
)

diagram <- ggplot() +
  geom_circle(
    data = circles,
    aes(x0 = x0, y0 = y0, r = r, fill = set),
    alpha = 0.25,
    color = "grey30",
    linewidth = 1
  ) +
  annotate("text", x = 0, y = 1.65, label = "MILK = CHILD", fontface = "bold", size = 6) +
  annotate("text", x = 0, y = -0.05, label = "LIFELINES", fontface = "bold", size = 6) +
  annotate("text", x = 0, y = 0.65, label = "16 shared phenotypes\n6 not in LIFELINES Cohort", size = 4) +
  coord_equal(xlim = c(-2.8, 2.8), ylim = c(-2.8, 2.8), expand = FALSE) +
  scale_fill_manual(values = c("A = B (22)" = "skyblue", "C (16)" = "salmon")) +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "Shared HMO Phenotypes Across Cohorts")
ggsave("results/euler_diag.png", plot = diagram,
       width = 14,
       height = 8,
       dpi = 300)

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

### metal summary tables
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
    filter(
      !is.na(P), P > 0, P < 5e-8,
      !is.na(CHR), !is.na(BP),
      !is.na(Direction),
      stringr::str_count(Direction, "\\?") <= 0 # 1 for 3-cohort run, 0 for milk/child run
    )
  
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
      Effect_Allele = Allele1,
      Other_Allele = Allele2,
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

# overall 3-cohort run
summ_df <- map_dfr(metal_files, summary_tbl)
summ_df$Phenotype <- sub("1\\.tbl$", "", summ_df$Phenotype)

library(gt)
summ_df %>%
  gt() %>%
  gtsave(filename = "results/lead_snps.html")

# secondary milk+child run
summ_df_mc <- map_dfr(milkchild_files, summary_tbl)
summ_df_mc$Phenotype <- sub("_child_milk1\\.tbl$", "", summ_df_mc$Phenotype)
summ_df_mc <- summ_df_mc %>%
  filter(Phenotype %in% c("DFLNH", "DFLNT", "DFLac", "DSLNH", "FDSLNH", "FLNH"))

summ_df_mc %>%
  gt() %>%
  gtsave(filename = "results/lead_snps_mc.html")

### forest plots
# first ensure the strands are not flipped
harmonize_beta <- function(beta, ea, oa, ref_ea, ref_oa) {
  beta   <- as.numeric(beta)
  ea     <- toupper(trimws(as.character(ea)))
  oa     <- toupper(trimws(as.character(oa)))
  ref_ea <- toupper(trimws(as.character(ref_ea)))
  ref_oa <- toupper(trimws(as.character(ref_oa)))
  
  comp <- function(x) chartr("ACGT", "TGCA", x)
  
  out <- list(beta_plot = NA_real_, match_type = "no_match")
  
  if (anyNA(c(beta, ea, oa, ref_ea, ref_oa))) {
    out$match_type <- "missing"
    return(out)
  }
  
  if (isTRUE(ea == ref_ea) && isTRUE(oa == ref_oa)) {
    out$beta_plot <- beta
    out$match_type <- "same"
    return(out)
  }
  
  if (isTRUE(ea == ref_oa) && isTRUE(oa == ref_ea)) {
    out$beta_plot <- -beta
    out$match_type <- "reversed"
    return(out)
  }
  
  ea_c <- comp(ea)
  oa_c <- comp(oa)
  
  if (isTRUE(ea_c == ref_ea) && isTRUE(oa_c == ref_oa)) {
    out$beta_plot <- beta
    out$match_type <- "complement"
    return(out)
  }
  
  if (isTRUE(ea_c == ref_oa) && isTRUE(oa_c == ref_ea)) {
    out$beta_plot <- -beta
    out$match_type <- "complement_reversed"
    return(out)
  }
  
  out
}

make_study_row <- function(df, id_col, beta_col, se_col, ea_col, oa_col,
                           obs_id, ref_ea, ref_oa, phenotype, cohort) {
  x <- df %>% filter(.data[[id_col]] == obs_id)
  
  if (nrow(x) == 0) return(NULL)
  
  x <- x %>% slice(1)
  
  hm <- harmonize_beta(
    beta   = x[[beta_col]][1],
    ea     = x[[ea_col]][1],
    oa     = x[[oa_col]][1],
    ref_ea = ref_ea,
    ref_oa = ref_oa
  )
  
  tibble(
    Phenotype = phenotype,
    SNP = x[[id_col]][1],
    Beta = hm$beta_plot,
    SE = x[[se_col]][1],
    Effect_Allele = x[[ea_col]][1],
    Other_Allele = x[[oa_col]][1],
    Match = hm$match_type,
    Cohort = cohort
  )
}

# now plot
forest_plt_df <- function(file_name, obs_id, ref_ea, ref_oa) {
  
  file1 <- paste0("data/milk/milk_", file_name, "_clean.mlma")
  file2 <- paste0("data/child/child_", file_name, "_hg38_cleaned.tsv")
  file3 <- paste0("data/lifeline/lifeline.", file_name, ".hg38.cleaned.tsv")
  
  study1 <- make_study_row(
    df = fread(file1),
    id_col = "SNP",
    beta_col = "b",
    se_col = "se",
    ea_col = "A1",
    oa_col = "A2",
    obs_id = obs_id,
    ref_ea = ref_ea,
    ref_oa = ref_oa,
    phenotype = file_name,
    cohort = "MILK"
  )
  
  study2 <- make_study_row(
    df = fread(file2),
    id_col = "rs_id",
    beta_col = "beta",
    se_col = "standard_error",
    ea_col = "effect_allele",
    oa_col = "other_allele",
    obs_id = obs_id,
    ref_ea = ref_ea,
    ref_oa = ref_oa,
    phenotype = file_name,
    cohort = "CHILD"
  )
  
  study3 <- make_study_row(
    df = fread(file3),
    id_col = "ID",
    beta_col = "BETA",
    se_col = "SE",
    ea_col = "ALT",
    oa_col = "REF",
    obs_id = obs_id,
    ref_ea = ref_ea,
    ref_oa = ref_oa,
    phenotype = file_name,
    cohort = "LIFELINE"
  )
  
  bind_rows(study1, study2, study3)
}

meta_ref <- summ_df %>%
  transmute(
    Phenotype,
    SNP,
    ref_ea = Effect_Allele,
    ref_oa = Other_Allele
  )

plot_forest_df <- pmap_dfr(
  list(meta_ref$Phenotype, meta_ref$SNP, meta_ref$ref_ea, meta_ref$ref_oa),
  forest_plt_df
)

summ_df$Cohort <- "Meta"
plot_forest_df <- bind_rows(
  plot_forest_df,
  summ_df %>%
    transmute(
      Phenotype,
      SNP,
      Beta,
      SE,
      Effect_Allele,
      Other_Allele,
      Cohort
    )
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


fp <- ggplot(plot_forest_df, aes(x = Beta, y = Cohort, color = Cohort)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
  geom_point(size = 2) +
  facet_wrap(~ Phenotype + SNP, scales = "free_y") +
  labs(x = "Estimate", y = "Cohort", title = "Beta Effect Sizes for Lead SNP Variants from Meta-analysis of All Cohorts") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/all_cohorts/forest_plot.png", plot = fp,
  width = 14,
  height = 8,
  dpi = 300)

### forest plots for child+milk cohorts
forest_plt_mc <- function(file_name, obs_id, ref_ea, ref_oa) {
  
  file1 <- paste0("data/milk/milk_", file_name, "_clean.mlma")
  file2 <- paste0("data/child/child_", file_name, "_hg38_cleaned.tsv")
  
  study1 <- make_study_row(
    df = fread(file1),
    id_col = "SNP",
    beta_col = "b",
    se_col = "se",
    ea_col = "A1",
    oa_col = "A2",
    obs_id = obs_id,
    ref_ea = ref_ea,
    ref_oa = ref_oa,
    phenotype = file_name,
    cohort = "MILK"
  )
  
  study2 <- make_study_row(
    df = fread(file2),
    id_col = "rs_id",
    beta_col = "beta",
    se_col = "standard_error",
    ea_col = "effect_allele",
    oa_col = "other_allele",
    obs_id = obs_id,
    ref_ea = ref_ea,
    ref_oa = ref_oa,
    phenotype = file_name,
    cohort = "CHILD"
  )
    
  bind_rows(study1, study2)
}

meta_ref_mc <- summ_df_mc %>%
  transmute(
    Phenotype,
    SNP,
    ref_ea = Effect_Allele,
    ref_oa = Other_Allele
  )

plot_forest_mc <- pmap_dfr(
  list(meta_ref_mc$Phenotype, meta_ref_mc$SNP, meta_ref_mc$ref_ea, meta_ref_mc$ref_oa),
  forest_plt_mc
)

summ_df_mc$Cohort <- "Meta"
plot_forest_mc <- bind_rows(
  plot_forest_mc,
  summ_df_mc %>%
    transmute(
      Phenotype,
      SNP,
      Beta,
      SE,
      Effect_Allele,
      Other_Allele,
      Cohort
    )
)

plot_forest_mc <- plot_forest_mc %>%
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

fp_mc <- ggplot(plot_forest_mc, aes(x = Beta, y = Cohort, color = Cohort)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
  geom_point(size = 2) +
  facet_wrap(~ Phenotype + SNP, scales = "free_y") +
  labs(x = "Estimate", y = "Cohort", title = "Beta Effect Sizes for Lead SNP Variants from Meta-analysis of Milk & Child Cohorts") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/milk_child_cohorts/forest_plot_mc.png", plot = fp_mc,
       width = 14,
       height = 8,
       dpi = 300)

### overlayed manhattan plot
read_tbl <- function(f) {
  dt <- fread(f)
  
  dt <- dt %>%
    separate(MarkerName,
             into = c("CHR", "BP"),
             sep = "_",
             extra = "drop",
             remove = FALSE) %>%
    mutate(
      trait = sub("^metal_(.+?)(?:_child_milk)?1\\.tbl$", "\\1", basename(f)),
      CHR = as.character(`CHR`),
      BP  = as.numeric(BP),
      P   = as.numeric(`P-value`)
    ) %>%
    filter(!is.na(P), P > 0, !is.na(CHR))
  
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
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10((5e-8)/6), linetype = "dashed", color = "blue") +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR) +
  labs(x = "Chromosome", y = expression(-log[10](P)), color = "Trait") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 3))
ggsave("results/overlay_manhattan_metal.png", plot = p,
       width = 12,
       height = 8,
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



