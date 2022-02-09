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
dir.plots <-  "plots/manhattan_plot/"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

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
sig.dmr.males <- readxl::read_xlsx(file.path("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/Michale/male_dist750/cnew.regions-p.bed_annotated.xlsx"))
sig.cpgs.in.dmr.males <- sig.dmr.males$cpgs_in_region %>% stringr::str_split(",") %>% unlist %>% unique()
length(sig.cpgs.in.dmr.males)

sig.dmr.females <- readxl::read_xlsx(file.path("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/Michale/female_dist750/cnew.regions-p.bed_annotated.xlsx"))
sig.cpgs.in.dmr.females <- sig.dmr.females$cpgs_in_region %>% stringr::str_split(",") %>% unlist %>% unique()
length(sig.cpgs.in.dmr.females)

cross.cpgs.males <- readr::read_csv(file.path("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/analysis_results/cross_meta_analysis/single_cpg/MALE_cross_tissue_meta_analysis_glm_using_AD_vs_CN_single_cpg_annotated.csv")
) %>% dplyr::filter(p < 10^-5  & valid_p == 6)

blood.meta_analysis.male <- readr::read_csv(
  "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/analysis_results/meta_analysis/Logistic_regression_model/AD_vs_CN/MALE_meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg_annotated.csv"
)
blood.meta_analysis.male <- blood.meta_analysis.male[,grep("Brain_sex_meta_analysis_Brain_sex_meta_analysis",colnames(blood.meta_analysis.male),invert = T)]

brain.meta_analysis.male <- readr::read_csv(
  "~/TBL Dropbox/Tiago Silva/AD_metaAnalysis_bySex/meta_analysis_single_cpg_results/single_cpg_male_meta_bacon_annot_df.csv"
)

table(cross.cpgs.males %in% blood.meta_analysis.male$cpg[blood.meta_analysis.male$pVal.final.bacon < 0.05])

sig.cpgs.in.meta.male <- intersect(brain.meta_analysis.male$cpg[brain.meta_analysis.male$pVal.final < 0.05],blood.meta_analysis.male$cpg[blood.meta_analysis.male$pVal.final.bacon < 0.05])
table(cross.cpgs.males$cpg %in% sig.cpgs.in.meta.male)

sig.cross.cpgs.males <- intersect(cross.cpgs.males$cpg,sig.cpgs.in.meta.male)
length(sig.cross.cpgs.males)


cross.cpgs.females <- readr::read_csv(file.path("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/analysis_results/cross_meta_analysis/single_cpg/FEMALE_cross_tissue_meta_analysis_glm_using_AD_vs_CN_single_cpg_with_NA_in_blood_or_brain_annotated.csv")
) %>% dplyr::filter(p < 10^-5 & valid_p == 6)


blood.meta_analysis.female <- readr::read_csv(
  "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/analysis_results/meta_analysis/Logistic_regression_model/AD_vs_CN/FEMALE_meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg_annotated.csv"
)
blood.meta_analysis.female <- blood.meta_analysis.female[,grep("Brain_sex_meta_analysis_",colnames(blood.meta_analysis.female),invert = T)]

brain.meta_analysis.female <- readr::read_csv(
  "~/TBL Dropbox/Tiago Silva/AD_metaAnalysis_bySex/meta_analysis_single_cpg_results/single_cpg_female_meta_bacon_annot_df.csv"
)

sig.cpgs.in.meta.female <- intersect(brain.meta_analysis.female$cpg[brain.meta_analysis.female$pVal.final < 0.05],blood.meta_analysis.female$cpg[blood.meta_analysis.female$pVal.final.bacon < 0.05])
table(cross.cpgs.females$cpg %in% sig.cpgs.in.meta.female)
sig.cross.cpgs.females <- intersect(cross.cpgs.females$cpg,sig.cpgs.in.meta.female)
length(sig.cross.cpgs.females)

dir.base <- "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex"
dir.results <- file.path(dir.base,"analysis_results")
dir.result.meta.analysis <- file.path(dir.results, "meta_analysis/Logistic_regression_model/")
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

sig.cpgs.in.dmr.females <- intersect(sig.cpgs.in.dmr.females,meta_df_AD_vs_CN.annotated.female$cpg)
sig.cpgs.in.dmr.males <- intersect(sig.cpgs.in.dmr.males,meta_df_AD_vs_CN.annotated.male$cpg)

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
plot <- ggmiami2(
  data = blood_male_female, 
  split_by = "Source", 
  split_at = "Blood male", 
  p = "pVal.final", 
  suggestive_line_color = "red",
  upper_ylab = "Blood male",
  lower_ylab = "Blood female",
  genome_line = NULL,
  top_n_hits = 10,
  suggestive_line_bottom = 10^-5,
  suggestive_line_upper = 10^-5,
  upper_labels_df = blood_labels_male,
  lower_labels_df = blood_labels_female
)


library(grid)
library(gtable)
ggplot2::ggsave(
  plot = gtable_add_padding(plot, unit(1, "cm")),
  filename = file.path(dir.plots,"blood_male_female_manhattan_plot.png"),
  width = 7,
  height = 6
)





