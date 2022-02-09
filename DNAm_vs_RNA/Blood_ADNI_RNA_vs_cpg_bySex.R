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
# Date: 21 July 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Article section: 
# Correlations between methylation levels of significant CpGs and DMRs in AD 
# with expressions of nearby genes
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Libs
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
library(dplyr)
library(SummarizedExperiment)
library(MethReg)


#-----------------------------------------------------------------------------
# Directories
#-----------------------------------------------------------------------------
dir.base <- "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/"
dir.results <- file.path(dir.base,"analysis_results")
dir.result.meta.analysis <- file.path(dir.results, "meta_analysis/Logistic_regression_model/")
dir.result.meta.analysis.ad <- file.path(dir.result.meta.analysis, "AD_vs_CN/")
dir.RNA_vs_DNAm <- "analysis_results/RNA_vs_DNAm/"
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)



#-------------------------------------------------------------------------------
# Data 
#-------------------------------------------------------------------------------
load("datasets/Aux/ADNI_matched_rna_dnam_residuals_bySex.rda")
all(rownames(metadata.dnam) == rownames(metadata.exp))
all(colnames(residuals.matched.exp) == colnames(residuals.matched.met))


#-------------------------------------------------------------------------------
# Aux functions
#-------------------------------------------------------------------------------
auxfunction <- function(
  row,
  exp = NULL,
  dnam = NULL,
  metadata.dnam = NULL
){
  
  rna.target <- exp[which(rownames(exp) == row$target), , drop = FALSE]
  met <- dnam[which(rownames(dnam) == as.character(row$regionID)), , drop = FALSE]
  
  data <- data.frame(
    "met.residual" = met %>% as.numeric(),
    "rna.target.residual" = rna.target %>% as.numeric(),
    "AD_Status" = factor(metadata.dnam$DX,levels = c("Dementia","CN"))
  )
  
  rlm <- MASS::rlm(
    rna.target.residual ~ met.residual + AD_Status, 
    data = data,
    psi = MASS::psi.bisquare,
    maxit = 100
  ) %>% summary %>% coef %>% data.frame
  
  degrees.freedom.value <- nrow(data) - 3
  rlm$pval <- 2 * (1 - pt( abs(rlm$t.value), df = degrees.freedom.value) )
  
  quant.pval <- rlm[-1,4,drop = FALSE] %>%
    t %>%
    as.data.frame()
  colnames(quant.pval) <- paste0("RLM_",colnames(quant.pval),"_pvalue")
  
  quant.estimate <- rlm[-1,1,drop = FALSE] %>%
    t %>%
    as.data.frame()
  colnames(quant.estimate) <- paste0("RLM_",colnames(quant.estimate),"_estimate")
  return(cbind(quant.pval, quant.estimate))
}


add_cpgs_from_and_do_fdr <- function(results, analysis = "FEMALE"){
  
  if(analysis == "FEMALE"){
    idx.ad.vs.cn <- which(results$probeID %in% sig.cpgs.females)
    idx.cross <- which(results$probeID %in% sig.cross.cpgs.females)
    idx.combp <- which(results$probeID %in% sig.cpgs.in.dmr.females)
    
    n.cpgs.ad <- length(sig.cpgs.females)
    n.cpgs.cross <- length(sig.cross.cpgs.females)
    n.cpgs.compb <- length(sig.cpgs.in.dmr.females)
    
    
  } else {
    idx.ad.vs.cn <- which(results$probeID %in% sig.cpgs.males)
    idx.cross <- which(results$probeID %in% sig.cross.cpgs.males)
    idx.combp <- which(results$probeID %in% sig.cpgs.in.dmr.males)
    
    n.cpgs.ad <- length(sig.cpgs.males)
    n.cpgs.cross <- length(sig.cross.cpgs.males)
    n.cpgs.compb <- length(sig.cpgs.in.dmr.males)
  }
  
  
  label.ad <- paste0(analysis," - ",n.cpgs.ad, " CpGs from AD vs CN meta-analysis")
  label.cross <- paste0(analysis," - ",n.cpgs.cross, " CpGs from cross-tissue meta-analysis")
  label.combp <- paste0(analysis," - ", n.cpgs.compb, " CpGs from comb-p meta-analysis")
  
  results <- results[c(idx.ad.vs.cn,idx.cross,idx.combp),]
  results$cpgs_from <- c(
    rep(label.ad, length(idx.ad.vs.cn)),
    rep(label.cross, length(idx.cross)),
    rep(label.combp, length(idx.combp))
  )
  
  
  results$RLM_met.residual_fdr <- c(
    p.adjust(results$RLM_met.residual_pvalue[which(results$cpgs_from %in% label.ad)],method = "fdr"),
    p.adjust(results$RLM_met.residual_pvalue[which(results$cpgs_from %in% label.cross)],method = "fdr"),
    p.adjust(results$RLM_met.residual_pvalue[which(results$cpgs_from %in% label.combp)],method = "fdr")
  )
  return(results)
}

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

sig.dmr.males <- readxl::read_xlsx(file.path("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/Michale/male_dist750/cnew.regions-p.bed_annotated.xlsx"))
sig.cpgs.in.dmr.males <- sig.dmr.males$cpgs_in_region %>% stringr::str_split(",") %>% unlist %>% unique()
length(sig.cpgs.in.dmr.males)

sig.dmr.females <- readxl::read_xlsx(file.path("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/Michale/female_dist750/cnew.regions-p.bed_annotated.xlsx"))
sig.cpgs.in.dmr.females <- sig.dmr.females$cpgs_in_region %>% stringr::str_split(",") %>% unlist %>% unique()
length(sig.cpgs.in.dmr.females)

#-----------------------------------------------------------------------------
# Select DMRs
#-----------------------------------------------------------------------------
cross.cpgs.males <- readr::read_csv(file.path("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/analysis_results/cross_meta_analysis/single_cpg/MALE_cross_tissue_meta_analysis_glm_using_AD_vs_CN_single_cpg_annotated.csv")
) %>% dplyr::filter(p < 10^-5  & valid_p == 6)

blood.meta_analysis.male <- readr::read_csv(
  "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/analysis_results/meta_analysis/Logistic_regression_model/AD_vs_CN/MALE_meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg_annotated.csv"
)
blood.meta_analysis.male <- blood.meta_analysis.male[,grep("Brain_sex_meta_analysis_Brain_sex_meta_analysis",colnames(blood.meta_analysis.male),invert = T)]

brain.meta_analysis.male <- readr::read_csv(
  "~/TBL Dropbox/Tiago Silva/AD_metaAnalysis_bySex/meta_analysis_single_cpg_results/single_cpg_male_meta_bacon_annot_df.csv"
)

intersect(brain.meta_analysis.male$cpg[brain.meta_analysis.male$pVal.final < 0.05],blood.meta_analysis.male$cpg[blood.meta_analysis.male$pVal.final.bacon < 0.05])
table(cross.cpgs.males %in% blood.meta_analysis.male$cpg[blood.meta_analysis.male$pVal.final.bacon < 0.05])

sig.cpgs.in.meta.male <- intersect(brain.meta_analysis.male$cpg[brain.meta_analysis.male$pVal.final < 0.05],blood.meta_analysis.male$cpg[blood.meta_analysis.male$pVal.final.bacon < 0.05])
table(cross.cpgs.males$cpg %in% sig.cpgs.in.meta.male)

sig.cross.cpgs.males <- intersect(cross.cpgs.males$cpg,sig.cpgs.in.meta.male)
length(sig.cross.cpgs.males)

cross.cpgs.females <- readr::read_csv(file.path("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/analysis_results/cross_meta_analysis/single_cpg/FEMALE_cross_tissue_meta_analysis_glm_using_AD_vs_CN_single_cpg_with_NA_in_blood_or_brain_annotated.csv")
) %>% dplyr::filter(p < 10^-5  & valid_p == 6)


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

sig.all.cpgs.male <- unique(c(sig.cpgs.in.dmr.males,sig.cpgs.males,sig.cross.cpgs.males))
sig.all.cpgs.female <- unique(c(sig.cpgs.in.dmr.females,sig.cpgs.females,sig.cross.cpgs.females))



#-------------------------------------------------------------------------------
# Promoter analysis, triplets using remap
#-------------------------------------------------------------------------------
dnam <- MethReg:::map_probes_to_regions(
  dnam = residuals.matched.met.male,
  genome = "hg19",
  arrayType = "EPIC",
  rm.masked.probes = FALSE
)
EPIC.hg19 <- MethReg:::get_met_probes_info(genome = "hg19", arrayType = "EPIC")
promoter.gene.dnam.pair <- MethReg::get_region_target_gene(
  rownames(dnam) %>% MethReg::make_granges_from_names(),
  method = "genes.promoter.overlap",
  genome = "hg19"
)

# 1 out of 5 cpgs in males are promoter
table((promoter.gene.dnam.pair$regionID %>% unique()) %in% make_names_from_granges(EPIC.hg19[sig.all.cpgs.male]))

nrow(promoter.gene.dnam.pair)
promoter.gene.dnam.pair <- promoter.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp.male))
nrow(promoter.gene.dnam.pair)
promoter.gene.dnam.pair$probeID <- names(EPIC.hg19)[match(promoter.gene.dnam.pair$regionID,make_names_from_granges(EPIC.hg19))]
promoter.gene.dnam.pair <- promoter.gene.dnam.pair %>% dplyr::filter(.data$regionID %in% rownames(dnam))
nrow(promoter.gene.dnam.pair)
results.promoter.analysis <- plyr::adply(
  .data = promoter.gene.dnam.pair,
  .margins = 1,
  .fun = function(row){
    auxfunction(row = row,exp = residuals.matched.exp.male,metadata.dnam = metadata.dnam.male,dnam = dnam)    
  }
)

results.promoter.analysis <- results.promoter.analysis %>% add_cpgs_from_and_do_fdr(analysis = "MALE")

results.promoter.analysis[results.promoter.analysis$RLM_met.residual_fdr < 0.05,]


#-------------------------------------------------------------------------------
# Distal analysis, triplets using remap
#-------------------------------------------------------------------------------
distal.gene.dnam.pair <- MethReg::get_region_target_gene(
  rownames(dnam) %>% MethReg::make_granges_from_names(),
  method = "nearby.genes",
  num.flanking.genes = 10,
  rm.promoter.regions.from.distal.linking = TRUE,
  genome = "hg19"
)

# 4 cpgs in are promoter
table((distal.gene.dnam.pair$regionID %>% unique()) %in% make_names_from_granges(EPIC.hg19[sig.all.cpgs.male]))

nrow(distal.gene.dnam.pair)
distal.gene.dnam.pair <- distal.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp.male))
nrow(distal.gene.dnam.pair)
distal.gene.dnam.pair$probeID <- names(EPIC.hg19)[match(distal.gene.dnam.pair$regionID,make_names_from_granges(EPIC.hg19))]
distal.gene.dnam.pair <- distal.gene.dnam.pair %>% dplyr::filter(.data$regionID %in% rownames(dnam))
nrow(distal.gene.dnam.pair)

results.distal.analysis <- plyr::adply(
  .data = distal.gene.dnam.pair,
  .margins = 1,
  .fun = function(row){
    auxfunction(row = row, exp = residuals.matched.exp.male, metadata.dnam = metadata.dnam.male, dnam = dnam)    
  }
)


results.distal.analysis <- results.distal.analysis %>% add_cpgs_from_and_do_fdr(analysis = "MALE")

results.distal.analysis[results.distal.analysis$RLM_met.residual_fdr < 0.05,]

results.promoter.analysis.male <- results.promoter.analysis
results.distal.analysis.male <- results.distal.analysis





#-------------------------------------------------------------------------------
# FEMALE 
#-------------------------------------------------------------------------------
meta_df_AD_vs_CN.annotated.sig.female <- readr::read_csv(
  file = paste0(
    dir.result.meta.analysis.ad,
    "FEMALE_meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg_annotated_pvalue_cut_off_0_05.csv"
  ),
  col_types = readr::cols()
)  %>% dplyr::filter(pVal.final.bacon < 1E-5) 



#-------------------------------------------------------------------------------
# Promoter analysis, triplets using remap
#-------------------------------------------------------------------------------
dnam <- MethReg:::map_probes_to_regions(
  dnam = residuals.matched.met.female,
  genome = "hg19",
  arrayType = "EPIC",
  rm.masked.probes = FALSE
)
EPIC.hg19 <- MethReg:::get_met_probes_info(genome = "hg19", arrayType = "EPIC")
promoter.gene.dnam.pair <- MethReg::get_region_target_gene(
  rownames(dnam) %>% MethReg::make_granges_from_names(),
  method = "genes.promoter.overlap",
  genome = "hg19"
)

# 2 out of in females are promoter
table((promoter.gene.dnam.pair$regionID %>% unique()) %in% make_names_from_granges(EPIC.hg19[sig.all.cpgs.female]))

nrow(promoter.gene.dnam.pair)
promoter.gene.dnam.pair <- promoter.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp.female))
nrow(promoter.gene.dnam.pair)
promoter.gene.dnam.pair$probeID <- names(EPIC.hg19)[match(promoter.gene.dnam.pair$regionID,make_names_from_granges(EPIC.hg19))]
promoter.gene.dnam.pair <- promoter.gene.dnam.pair %>% dplyr::filter(.data$regionID %in% rownames(dnam))
nrow(promoter.gene.dnam.pair)
results.promoter.analysis <- plyr::adply(
  .data = promoter.gene.dnam.pair,
  .margins = 1,
  .fun = function(row){
    auxfunction(row = row,exp = residuals.matched.exp.female,metadata.dnam = metadata.dnam.female,dnam = dnam)    
  }
)

results.promoter.analysis <- results.promoter.analysis %>% add_cpgs_from_and_do_fdr(analysis = "FEMALE")

results.promoter.analysis[results.promoter.analysis$RLM_met.residual_fdr < 0.05,]


#-------------------------------------------------------------------------------
# Distal analysis, triplets using remap
#-------------------------------------------------------------------------------
distal.gene.dnam.pair <- MethReg::get_region_target_gene(
  rownames(dnam) %>% MethReg::make_granges_from_names(),
  method = "nearby.genes",
  num.flanking.genes = 10,
  rm.promoter.regions.from.distal.linking = TRUE,
  genome = "hg19"
)

# 18 cpgs in are promoter
table((distal.gene.dnam.pair$regionID %>% unique()) %in% make_names_from_granges(EPIC.hg19[sig.all.cpgs.female]))

nrow(distal.gene.dnam.pair)
distal.gene.dnam.pair <- distal.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp.female))
nrow(distal.gene.dnam.pair)
distal.gene.dnam.pair$probeID <- names(EPIC.hg19)[match(distal.gene.dnam.pair$regionID,make_names_from_granges(EPIC.hg19))]
distal.gene.dnam.pair <- distal.gene.dnam.pair %>% dplyr::filter(.data$regionID %in% rownames(dnam))
nrow(distal.gene.dnam.pair)

results.distal.analysis <- plyr::adply(
  .data = distal.gene.dnam.pair,
  .margins = 1,
  .fun = function(row){
    auxfunction(row = row, exp = residuals.matched.exp.female, metadata.dnam = metadata.dnam.female, dnam = dnam)    
  }
)


results.distal.analysis <- results.distal.analysis %>% add_cpgs_from_and_do_fdr("FEMALE")

results.distal.analysis[results.distal.analysis$RLM_met.residual_fdr < 0.05,]

results.promoter.analysis.female <- results.promoter.analysis
results.distal.analysis.female <- results.distal.analysis



#-------------------------------------------------------------------------------
# Save results
#-------------------------------------------------------------------------------

writexl::write_xlsx(
  list(
    "MALE_Promoter" = results.promoter.analysis.male,
    "MALE_Distal_10_up_10_down" = results.distal.analysis.male,
    "FEMALE_Promoter" = results.promoter.analysis.female,
    "FEMALE_Distal_10_up_10_down" = results.distal.analysis.female
  ),
  path = file.path(dir.RNA_vs_DNAm,"FEMALE_MALE_Blood_ADNI_cpg_Target_vs_DNAm.xlsx")
)
