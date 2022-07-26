#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Article:
# Integrative meta-analysis of epigenome-wide association studies
# identifies genomic and
# epigenomics differences in the brain and the blood in Alzheimer’s disease
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Authors: 
# - Tiago C. silva
# - Lily Wang
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Date: 12 July 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.plots <-  "~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/plots/manhattan_plot/withSmoke"
dir.create(dir.plots,recursive = TRUE,showWarnings = FALSE)

#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------
# install_github("juliedwhite/miamiplot", build_vignettes = TRUE)
library(miamiplot)
library(dplyr)
library(dplyr)
source("code/manhattan_plot/ggmiami.R")

#-------------------------------------------------------------------------------
# Get brain meta-analysis results data
#-------------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Select cpgs
#-----------------------------------------------------------------------------
# FEMALE CpGs with P<1E- 5 in AD vs. CN comparison
sig.cpgs.females <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/DRAFT-FIGURES-TABLES_1-11-2021/Supp Table 2 sig_female_male_cpgs_AD_vs_CN.xlsx",skip = 3,n_max = 25)
sig.cpgs.females <- sig.cpgs.females[-1,]
header <- colnames(sig.cpgs.females)
sig.cpgs.females <- sig.cpgs.females$CpG
length(sig.cpgs.females)


# MALE CpGs with P<1E- 5 in AD vs. CN comparison
sig.cpgs.males <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/DRAFT-FIGURES-TABLES_1-11-2021/Supp Table 2 sig_female_male_cpgs_AD_vs_CN.xlsx", skip = 30)
colnames(sig.cpgs.males) <- header
sig.cpgs.males <- sig.cpgs.males$CpG
length(sig.cpgs.males)

#-----------------------------------------------------------------------------
# Select DMRs
#-----------------------------------------------------------------------------


# ----------------------------------------------------------------------------
# Directly load sig cpgs from list
# ----------------------------------------------------------------------------

list.sig.all.cpgs.female <- readxl::read_xlsx(
  "analysis_results/others/sig_cpg_table/FEMALE_all_significant_cpgs_updated.xlsx"
)
sig.all.cpgs.female <- list.sig.all.cpgs.female$cpgs
sig.cpgs.females <- list.sig.all.cpgs.female[
  lapply(str_split(list.sig.all.cpgs.female$Sources, ","), function(l) any(l %in% "AD vs CN")) %>% unlist(),
] %>% pull(cpgs)
sig.cpgs.in.dmr.females <- list.sig.all.cpgs.female[
  lapply(str_split(list.sig.all.cpgs.female$Sources, ","), function(l) any(l %in% "comb-p")) %>% unlist(),
] %>% pull(cpgs)
sig.cross.cpgs.females <- list.sig.all.cpgs.female[
  lapply(str_split(list.sig.all.cpgs.female$Sources, ","), function(l) any(l %in% "cross-tissue")) %>% unlist(),
] %>% pull(cpgs)

list.sig.all.cpgs.male <- readxl::read_xlsx(
  "analysis_results/others/sig_cpg_table/MALE_all_significant_cpgs.xlsx"
)
sig.all.cpgs.male <- list.sig.all.cpgs.male$cpgs
sig.cpgs.males <- list.sig.all.cpgs.male %>% filter(
  Sources %in% "AD vs CN"
) %>% pull(cpgs)
sig.cpgs.in.dmr.males <- list.sig.all.cpgs.male %>% filter(
  Sources %in% "comb-p"
) %>% pull(cpgs)
sig.cross.cpgs.males <- list.sig.all.cpgs.male %>% filter(
  Sources %in% "cross-tissue"
) %>% pull(cpgs)
# ----------------------------------------------------------------------------


dir.base <- "~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex"
dir.results <- file.path(dir.base,"analysis_results")
dir.result.meta.analysis <- file.path(dir.results, "meta_analysis/Logistic_regression_model/withSmokePrbs")
dir.result.meta.analysis.ad <- file.path(dir.result.meta.analysis, "AD_vs_CN/")

# data for comb-p
meta_df_AD_vs_CN.annotated.male <- readr::read_csv(
  file = paste0(
    dir.result.meta.analysis.ad,
    "MALE_meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg_annotated.csv"
  )
) %>% as.data.frame()

meta_df_AD_vs_CN.annotated.female <- readr::read_csv(
  file = paste0(
    dir.result.meta.analysis.ad,
    "FEMALE_meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg_annotated.csv"
  )
) %>% as.data.frame()


cpgs.blood.female <- data.frame(
  "CpG" = meta_df_AD_vs_CN.annotated.female$cpg,
  "Source" = "Blood female",
  "pVal.final" = meta_df_AD_vs_CN.annotated.female$pVal.final.bacon,
  "pos" = meta_df_AD_vs_CN.annotated.female$pos,
  "chr" = meta_df_AD_vs_CN.annotated.female$chr,
  "UCSC_RefGene_Name" = meta_df_AD_vs_CN.annotated.female$UCSC_RefGene_Name
)
cpgs.blood.female$chr <- as.numeric(gsub("chr","",cpgs.blood.female$chr))
cpgs.blood.female$chr <- ifelse(
  is.na(cpgs.blood.female$chr), 23, cpgs.blood.female$chr
)

cpgs.blood.male <- data.frame(
  "CpG" = meta_df_AD_vs_CN.annotated.male$cpg,
  "Source" = "Blood male",
  "pVal.final" = meta_df_AD_vs_CN.annotated.male$pVal.final.bacon,
  "pos" = meta_df_AD_vs_CN.annotated.male$pos,
  "chr" = meta_df_AD_vs_CN.annotated.male$chr,
  "UCSC_RefGene_Name" = meta_df_AD_vs_CN.annotated.male$UCSC_RefGene_Name
)
cpgs.blood.male$chr <- as.numeric(gsub("chr","",cpgs.blood.male$chr))
cpgs.blood.male$chr <- ifelse(
  is.na(cpgs.blood.male$chr), 23, cpgs.blood.male$chr
)

blood_male_female <- rbind(cpgs.blood.male,cpgs.blood.female)
#-------------------------------------------------------------------------------
# Get brain meta-analysis results data
#-------------------------------------------------------------------------------
library(miamiplot)
plot_data <- prep_miami_data(
  data = blood_male_female, split_by = "Source", split_at = "Blood male", p = "pVal.final"
)

#-------------------------------------------------------------------------------
# Label
#-------------------------------------------------------------------------------
blood_labels_female <- plot_data$lower %>%
  filter(CpG %in% sig.cpgs.females) %>%
  select(rel_pos, logged_p, UCSC_RefGene_Name) %>%
  dplyr::rename(label = UCSC_RefGene_Name) %>%
  tidyr::separate_rows(label, sep = ";") %>%
  unique() %>%
  arrange(desc(logged_p)) %>% na.omit()
blood_labels_female <- blood_labels_female[blood_labels_female$label != "",]
blood_labels_female <- blood_labels_female[!duplicated(blood_labels_female$label),]

blood_labels_male <- plot_data$upper %>%
  filter(CpG %in% sig.cpgs.males) %>%
  select(rel_pos, logged_p, UCSC_RefGene_Name) %>%
  dplyr::rename(label = UCSC_RefGene_Name) %>%
  tidyr::separate_rows(label, sep = ";") %>%
  unique() %>%
  arrange(desc(logged_p)) %>% na.omit()
blood_labels_male <- blood_labels_male[blood_labels_male$label != "",]
blood_labels_male <- blood_labels_male[!duplicated(blood_labels_male$label),]


#-------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------
plot <- ggmiami(
  data = blood_male_female, 
  split_by = "Source", 
  split_at = "Blood male", 
  p = "pVal.final", 
  suggestive_line_color = "red",
  upper_ylab = "Blood male",
  lower_ylab = "Blood female",
  genome_line = NULL,
  top_n_hits = 10,
  suggestive_line = 10^-5,
  upper_labels_df = blood_labels_male,
  lower_labels_df = blood_labels_female
)


library(grid)
library(gtable)
ggplot2::ggsave(
  plot = gtable_add_padding(plot, unit(1, "cm")),
  filename = file.path(dir.plots,"blood_male_female_manhattan_plot.pdf"),
  width = 7,
  height = 6
)





