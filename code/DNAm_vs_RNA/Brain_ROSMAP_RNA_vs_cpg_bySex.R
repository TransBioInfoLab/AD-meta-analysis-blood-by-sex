#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Article:
# Integrative meta-analysis of epigenome-wide association studies
# identifies genomic and
# epigenomics differences in the brain and the blood in Alzheimer’s disease
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Authors: 
# - Tiago C. silva
# - Wei Zhang
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
library(coMethDMR)
#-----------------------------------------------------------------------------
# Analysis: target gene ~ CpG  uisng ROSMAP data
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Select cpgs
#-----------------------------------------------------------------------------
# FEMALE CpGs with P<1E- 5 in AD vs. CN comparison
sig.cpgs.females <- readxl::read_xlsx("~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/DRAFT-FIGURES-TABLES-REVISION_7-1-2022/Supp Table 2 sig_female_male_cpgs_AD_vs_CN_updated_final_LWedit - FINAL.xlsx",skip = 3,n_max = 24)
sig.cpgs.females <- sig.cpgs.females[-1,]
header <- colnames(sig.cpgs.females)
sig.cpgs.females <- sig.cpgs.females$cpg
length(sig.cpgs.females)


# MALE CpGs with P<1E- 5 in AD vs. CN comparison
sig.cpgs.males <- readxl::read_xlsx("~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/DRAFT-FIGURES-TABLES-REVISION_7-1-2022/Supp Table 2 sig_female_male_cpgs_AD_vs_CN_updated_final_LWedit - FINAL.xlsx", skip = 29)
colnames(sig.cpgs.males) <- header
sig.cpgs.males <- sig.cpgs.males$cpg
length(sig.cpgs.males)

sig.dmr.males <- readxl::read_xlsx(file.path("~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/Michale/male_750/cnew.regions-p.bed_annotated_with_directions_final_z_p_filtered.xlsx"))
sig.cpgs.in.dmr.males <- sig.dmr.males$cpgs_in_region %>% stringr::str_split(",") %>% unlist %>% unique()
length(sig.cpgs.in.dmr.males)

sig.dmr.females <- readxl::read_xlsx(file.path("~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/Michale/female_750/cnew.regions-p.bed_annotated_with_directions_final_z_p_filtered.xlsx"))
sig.cpgs.in.dmr.females <- sig.dmr.females$cpgs_in_region %>% stringr::str_split(",") %>% unlist %>% unique()
length(sig.cpgs.in.dmr.females)

#-----------------------------------------------------------------------------
# Select DMRs
#-----------------------------------------------------------------------------
cross.cpgs.males <- readr::read_csv(file.path("~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/analysis_results/cross_meta_analysis/withSmoke/MALE_cross_tissue_meta_analysis_glm_using_AD_vs_CN_single_cpg_annotated.csv")
) %>% dplyr::filter(p < 10^-5 & valid_p == 6)

blood.meta_analysis.male <- readr::read_csv(
  "~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/analysis_results/meta_analysis/Logistic_regression_model/withSmokePrbs/AD_vs_CN/MALE_meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg_annotated.csv"
)
blood.meta_analysis.male <- blood.meta_analysis.male[,grep("Brain_sex_meta_analysis_Brain_sex_meta_analysis",colnames(blood.meta_analysis.male),invert = T)]

brain.meta_analysis.male <- readr::read_csv(
  "~/TBL Dropbox/Wei Zhang/AD_metaAnalysis_bySex/meta_analysis_single_cpg_results/single_cpg_male_meta_bacon_annot_df.csv"
)

intersect(brain.meta_analysis.male$cpg[brain.meta_analysis.male$pVal.final < 0.05],blood.meta_analysis.male$cpg[blood.meta_analysis.male$pVal.final.bacon < 0.05])
table(cross.cpgs.males %in% blood.meta_analysis.male$cpg[blood.meta_analysis.male$pVal.final.bacon < 0.05])

sig.cpgs.in.meta.male <- intersect(brain.meta_analysis.male$cpg[brain.meta_analysis.male$pVal.final < 0.05],blood.meta_analysis.male$cpg[blood.meta_analysis.male$pVal.final.bacon < 0.05])
table(cross.cpgs.males$cpg %in% sig.cpgs.in.meta.male)

sig.cross.cpgs.males <- intersect(cross.cpgs.males$cpg,sig.cpgs.in.meta.male)
length(sig.cross.cpgs.males)

cross.cpgs.females <- readr::read_csv(file.path("~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/analysis_results/cross_meta_analysis/withSmoke/FEMALE_cross_tissue_meta_analysis_glm_using_AD_vs_CN_single_cpg_with_NA_in_blood_or_brain_annotated.csv")
) %>% dplyr::filter(p < 10^-5  & valid_p == 6)


blood.meta_analysis.female <- readr::read_csv(
  "~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/analysis_results/meta_analysis/Logistic_regression_model/withSmokePrbs/AD_vs_CN/FEMALE_meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg_annotated.csv"
)
blood.meta_analysis.female <- blood.meta_analysis.female[,grep("Brain_sex_meta_analysis_",colnames(blood.meta_analysis.female),invert = T)]

brain.meta_analysis.female <- readr::read_csv(
  "~/TBL Dropbox/Wei Zhang/AD_metaAnalysis_bySex/meta_analysis_single_cpg_results/single_cpg_female_meta_bacon_annot_df.csv"
)

sig.cpgs.in.meta.female <- intersect(brain.meta_analysis.female$cpg[brain.meta_analysis.female$pVal.final < 0.05],blood.meta_analysis.female$cpg[blood.meta_analysis.female$pVal.final.bacon < 0.05])
table(cross.cpgs.females$cpg %in% sig.cpgs.in.meta.female)
sig.cross.cpgs.females <- intersect(cross.cpgs.females$cpg,sig.cpgs.in.meta.female)
length(sig.cross.cpgs.females)

dir.base <- "~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex"
dir.results <- file.path(dir.base,"analysis_results")
dir.result.meta.analysis <- file.path(dir.results, "meta_analysis/Logistic_regression_model/withSmokePrbs/")
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

#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
# ROSMAP DNAm data
#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
# Data not included
load("coMethDMR_metaAnalysis/DNAm_RNA/data/matched_data.rda")
dim(matched.dnam)
dim(matched.exp)

matched.phenotype.female <- matched.phenotype[matched.phenotype$sex == "Female",]
matched.phenotype.male <- matched.phenotype[matched.phenotype$sex == "Male",]

matched.exp.male <- matched.exp[,match(matched.phenotype.male$rnaseq_id,gsub("_[0-9]*$","",colnames(matched.exp)))]
matched.exp.male <- matched.exp.male[rowSums(matched.exp.male) > 0,]

matched.exp.female <- matched.exp[,match(matched.phenotype.female$rnaseq_id,gsub("_[0-9]*$","",colnames(matched.exp)))]
matched.exp.female <- matched.exp.female[rowSums(matched.exp.female) > 0,]


matched.dnam.male <- matched.dnam[,match(matched.phenotype.male$Sample,colnames(matched.dnam))]
matched.dnam.female <- matched.dnam[,match(matched.phenotype.female$Sample,colnames(matched.dnam))]


# 1) remove confounding effects in DNAm data: 
resid_met.female <- coMethDMR:::GetResiduals(
  dnam = matched.dnam.female[rownames(matched.dnam.female) %in% sig.all.cpgs.female,],
  betaToM = TRUE, #converts to Mvalues for fitting linear model 
  pheno_df = matched.phenotype.female,
  covariates_char = c("Sample_Plate", "prop.neuron", "batch","age_death"), 
  nCores_int = 1,
  progressbar = TRUE  
)

resid_dnam.female<- MethReg:::map_probes_to_regions(
  dnam = resid_met.female,
  genome = "hg19",
  arrayType = "EPIC",
  rm.masked.probes = FALSE
)


resid_met.male <- coMethDMR:::GetResiduals(
  dnam = matched.dnam.male[rownames(matched.dnam.male) %in% sig.all.cpgs.male,],
  betaToM = TRUE, #converts to Mvalues for fitting linear model 
  pheno_df = matched.phenotype.male,
  covariates_char = c("Sample_Plate", "prop.neuron", "batch","age_death"), 
  nCores_int = 1,
  progressbar = TRUE  
)

resid_dnam.male <- MethReg:::map_probes_to_regions(
  dnam = resid_met.male,
  genome = "hg19",
  arrayType = "EPIC",
  rm.masked.probes = FALSE
)
#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
# Get linked genes
#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
EPIC.hg19 <- MethReg:::get_met_probes_info(genome = "hg19", arrayType = "EPIC")

promoter.gene.dnam.pair.female <- MethReg::get_region_target_gene(
  rownames(resid_dnam.female) %>% MethReg::make_granges_from_names(),
  method = "genes.promoter.overlap",
  genome = "hg19"
)

promoter.gene.dnam.pair.female$probeID <- names(EPIC.hg19)[match(promoter.gene.dnam.pair.female$regionID,make_names_from_granges(EPIC.hg19))]

distal.gene.dnam.pair.female <- MethReg::get_region_target_gene(
  rownames(resid_dnam.female) %>% MethReg::make_granges_from_names(),
  method = "nearby.genes",
  num.flanking.genes = 10,
  rm.promoter.regions.from.distal.linking = TRUE,
  genome = "hg19"
)

distal.gene.dnam.pair.female$probeID <- names(EPIC.hg19)[match(distal.gene.dnam.pair.female$regionID,make_names_from_granges(EPIC.hg19))]


window.gene.dnam.pair.famale <- MethReg::get_region_target_gene(
  rownames(resid_dnam.female) %>% MethReg::make_granges_from_names(),
  method = "window",
  genome = "hg19",
  window.size = 500 * 10^3,
  rm.promoter.regions.from.distal.linking = FALSE
)
window.gene.dnam.pair.famale$probeID <- names(EPIC.hg19)[match(window.gene.dnam.pair.famale$regionID,make_names_from_granges(EPIC.hg19))]

#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
# MALE: Get linked genes
#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=

promoter.gene.dnam.pair.male <- MethReg::get_region_target_gene(
  rownames(resid_dnam.male) %>% MethReg::make_granges_from_names(),
  method = "genes.promoter.overlap",
  genome = "hg19"
)

promoter.gene.dnam.pair.male$probeID <- names(EPIC.hg19)[match(promoter.gene.dnam.pair.male$regionID,make_names_from_granges(EPIC.hg19))]

distal.gene.dnam.pair.male <- MethReg::get_region_target_gene(
  rownames(resid_dnam.male) %>% MethReg::make_granges_from_names(),
  method = "nearby.genes",
  num.flanking.genes = 10,
  rm.promoter.regions.from.distal.linking = TRUE,
  genome = "hg19"
)

distal.gene.dnam.pair.male$probeID <- names(EPIC.hg19)[match(distal.gene.dnam.pair.male$regionID,make_names_from_granges(EPIC.hg19))]


window.gene.dnam.pair.male <- MethReg::get_region_target_gene(
  rownames(resid_dnam.male) %>% MethReg::make_granges_from_names(),
  method = "window",
  genome = "hg19",
  window.size = 500 * 10^3,
  rm.promoter.regions.from.distal.linking = FALSE
)
window.gene.dnam.pair.male$probeID <- names(EPIC.hg19)[match(window.gene.dnam.pair.male$regionID,make_names_from_granges(EPIC.hg19))]



#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
# ROSMAP Gene expression data
#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
matched.exp.log2.female <- log2(matched.exp.female + 1) # + 1 is required otherwise -INF
markers.female <-
  t(matched.exp.log2.female[c(
    "ENSG00000111674",
    "ENSG00000129226",
    "ENSG00000131095",
    "ENSG00000205927",
    "ENSG00000174059"
  ), ])
colnames(markers.female) <- c(
  "markers_ENO2",
  "markers_OLIG2",
  "markers_CD34",
  "markers_CD68",
  "markers_GFAP"
)

matched.exp.log2.female <- matched.exp.log2.female[ rownames(matched.exp.log2.female) %in% 
                                                      c(distal.gene.dnam.pair.female$target,
                                                        promoter.gene.dnam.pair.female$target,
                                                        window.gene.dnam.pair.famale$target),]

matched.phenotype.female$rnaseq_id  <- map$rnaseq_id[match(matched.phenotype.female$Sample,map$mwas_id)]

residuals.matched.exp.female <- plyr::adply(
  .data = matched.exp.log2.female,
  .margins = 1, 
  function(row){
    val <- t(row)
    colnames(val) <- "val"
    dat <- cbind(
      val, 
      matched.phenotype.female,
      markers.female
    )
    dat$val <- as.numeric(dat$val)
    fitE <- lm(
      "val ~ age_death + markers_ENO2 + markers_OLIG2 + markers_CD34 + markers_CD68 + markers_GFAP", 
      data = dat, 
      na.action = na.exclude
    )
    residuals(fitE)
  }, .progress = "time",
  .parallel = FALSE)
rownames(residuals.matched.exp.female) <- rownames(matched.exp.log2.female)


#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
# MALE ROSMAP Gene expression data
#-=-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-==-=-=-=-=--=-=
matched.exp.log2.male <- log2(matched.exp.male + 1) # + 1 is required otherwise -INF
markers.male <-
  t(matched.exp.log2.male[c(
    "ENSG00000111674",
    "ENSG00000129226",
    "ENSG00000131095",
    "ENSG00000205927",
    "ENSG00000174059"
  ), ])
colnames(markers.male) <- c(
  "markers_ENO2",
  "markers_OLIG2",
  "markers_CD34",
  "markers_CD68",
  "markers_GFAP"
)

matched.exp.log2.male <- matched.exp.log2.male[ rownames(matched.exp.log2.male) %in% 
                                                      c(distal.gene.dnam.pair.male$target,
                                                        promoter.gene.dnam.pair.male$target,
                                                        window.gene.dnam.pair.male$target),]

matched.phenotype.male$rnaseq_id  <- map$rnaseq_id[match(matched.phenotype.male$Sample,map$mwas_id)]

residuals.matched.exp.male <- plyr::adply(
  .data = matched.exp.log2.male,
  .margins = 1, 
  function(row){
    val <- t(row)
    colnames(val) <- "val"
    dat <- cbind(
      val, 
      matched.phenotype.male,
      markers.male
    )
    dat$val <- as.numeric(dat$val)
    fitE <- lm(
      "val ~ age_death + markers_ENO2 + markers_OLIG2 + markers_CD34 + markers_CD68 + markers_GFAP", 
      data = dat, 
      na.action = na.exclude
    )
    residuals(fitE)
  }, .progress = "time",
  .parallel = FALSE)
rownames(residuals.matched.exp.male) <- rownames(matched.exp.log2.male)


#-------------------------------------------------------------------------------
# Aux functions
#-------------------------------------------------------------------------------
auxfunction<- function(row, analysis = "FEMALE"){
  
  if(analysis == "FEMALE"){
    rna.target <- residuals.matched.exp.female[which(rownames(residuals.matched.exp.female) == row$target), , drop = FALSE]
    met <- resid_dnam.female[which(rownames(resid_dnam.female) == as.character(row$regionID)), , drop = FALSE]
    
    data <- data.frame(
      "met.residual" = met %>% as.numeric(),
      "rna.target.residual" = rna.target %>% as.numeric(),
      Braak_stage = matched.phenotype.female$braaksc %>% as.numeric()
    )
  } else {
    rna.target <- residuals.matched.exp.male[which(rownames(residuals.matched.exp.male) == row$target), , drop = FALSE]
    met <- resid_dnam.male[which(rownames(resid_dnam.male) == as.character(row$regionID)), , drop = FALSE]
    
    data <- data.frame(
      "met.residual" = met %>% as.numeric(),
      "rna.target.residual" = rna.target %>% as.numeric(),
      Braak_stage = matched.phenotype.male$braaksc %>% as.numeric()
    )
  }
  
  rlm <- MASS::rlm(
    rna.target.residual ~ met.residual + Braak_stage, 
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
   # idx.ad.vs.cn <- which(results$probeID %in% sig.cpgs.females)
    idx.cross <- which(results$probeID %in% sig.cross.cpgs.females)
   # idx.combp <- which(results$probeID %in% sig.cpgs.in.dmr.females)
    
  #  n.cpgs.ad <- length(sig.cpgs.females)
    n.cpgs.cross <- length(sig.cross.cpgs.females)
   # n.cpgs.compb <- length(sig.cpgs.in.dmr.females)
    
    
  } else {
    #idx.ad.vs.cn <- which(results$probeID %in% sig.cpgs.males)
    idx.cross <- which(results$probeID %in% sig.cross.cpgs.males)
    #idx.combp <- which(results$probeID %in% sig.cpgs.in.dmr.males)
    
    #n.cpgs.ad <- length(sig.cpgs.males)
    n.cpgs.cross <- length(sig.cross.cpgs.males)
    #n.cpgs.compb <- length(sig.cpgs.in.dmr.males)
  }
  
  
  #label.ad <- paste0(analysis," - ",n.cpgs.ad, " CpGs from AD vs CN meta-analysis")
  label.cross <- paste0(analysis," - ",n.cpgs.cross, " CpGs from cross-tissue meta-analysis")
  #label.combp <- paste0(analysis," - ", n.cpgs.compb, " CpGs from comb-p meta-analysis")
  
  results <- results[c(idx.cross),]
  results$cpgs_from <- c(
   # rep(label.ad, length(idx.ad.vs.cn)),
    rep(label.cross, length(idx.cross))
   # rep(label.combp, length(idx.combp))
  )
  
  
  results$RLM_met.residual_fdr <- c(
   # p.adjust(results$RLM_met.residual_pvalue[which(results$cpgs_from %in% label.ad)],method = "fdr"),
    p.adjust(results$RLM_met.residual_pvalue[which(results$cpgs_from %in% label.cross)],method = "fdr")
   # p.adjust(results$RLM_met.residual_pvalue[which(results$cpgs_from %in% label.combp)],method = "fdr")
  )
  return(results)
}


#-------------------------------------------------------------------------------
# Promoter analysis, triplets using remap
#-------------------------------------------------------------------------------

promoter.gene.dnam.pair.female <- promoter.gene.dnam.pair.female %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp.female))
results.promoter.analysis.female <- plyr::adply(promoter.gene.dnam.pair.female,.margins = 1,.fun = function(row) {auxfunction(row = row, analysis = "FEMALE")})

# Where did the cpg come from ? 

results.promoter.analysis.female <- results.promoter.analysis.female %>% add_cpgs_from_and_do_fdr(analysis = "FEMALE")

results.promoter.analysis.female[results.promoter.analysis.female$RLM_met.residual_fdr < 0.05,]

#-------------------------------------------------------------------------------
# Distal analysis, triplets using remap
#-------------------------------------------------------------------------------
nrow(distal.gene.dnam.pair.female) #  1068
distal.gene.dnam.pair.female <- distal.gene.dnam.pair.female %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp.female))
nrow(distal.gene.dnam.pair.female) # 863

results.distal.analysis.female <- plyr::adply(distal.gene.dnam.pair.female,.margins = 1,.fun = function(row) {auxfunction(row = row, analysis = "FEMALE")})

results.distal.analysis.female <- results.distal.analysis.female %>% add_cpgs_from_and_do_fdr(analysis = "FEMALE")

results.distal.analysis.female[results.distal.analysis.female$RLM_met.residual_fdr < 0.05,]


#-------------------------------------------------------------------------------
# Promoter analysis, triplets using remap
#-------------------------------------------------------------------------------

promoter.gene.dnam.pair.male <- promoter.gene.dnam.pair.male %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp.male))
results.promoter.analysis.male <- plyr::adply(promoter.gene.dnam.pair.male,.margins = 1,.fun = function(row) {auxfunction(row = row, analysis = "MALE")})
# Where did the cpg come from ? 

results.promoter.analysis.male <- results.promoter.analysis.male %>% add_cpgs_from_and_do_fdr(analysis = "MALE")

results.promoter.analysis.male[results.promoter.analysis.male$RLM_met.residual_fdr < 0.05,]

#-------------------------------------------------------------------------------
# Distal analysis, triplets using remap
#-------------------------------------------------------------------------------
nrow(distal.gene.dnam.pair.male) #  1040
distal.gene.dnam.pair.male <- distal.gene.dnam.pair.male %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp.male))
nrow(distal.gene.dnam.pair.male) # 753

results.distal.analysis.male <- plyr::adply(distal.gene.dnam.pair.male,.margins = 1,.fun = function(row) {auxfunction(row = row, analysis = "MALE")})

results.distal.analysis.male <- results.distal.analysis.male %>% add_cpgs_from_and_do_fdr(analysis = "MALE")

results.distal.analysis.male[results.distal.analysis.male$RLM_met.residual_fdr < 0.05,]


#-------------------------------------------------------------------------------
# Save resuls
#-------------------------------------------------------------------------------
path.RNA_vs_DNAm <- "~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex//analysis_results/RNA_vs_DNAm/withSmoke"
writexl::write_xlsx(
  list(
    "MALE Promoter" = results.promoter.analysis.male,
    "MALE Distal_10_up_10_down" = results.distal.analysis.male,
    "FEMALE Promoter" = results.promoter.analysis.female,
    "FEMALE Distal_10_up_10_down" = results.distal.analysis.female
    
  ),
  path = file.path(path.RNA_vs_DNAm,"FEMALE_MALE_Brain_cpg_Target_vs_DNAm.xlsx")
)


