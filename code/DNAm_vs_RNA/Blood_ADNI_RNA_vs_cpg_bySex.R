#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Article:
# Distinct sex-specific DNA methylation differences in Alzheimer’s disease
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
library(MethReg)


#-----------------------------------------------------------------------------
# Directories
#-----------------------------------------------------------------------------
dir.base <- "~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/"
dir.results <- file.path(dir.base,"analysis_results")
dir.result.meta.analysis <- file.path(dir.results, "meta_analysis/Logistic_regression_model/withSmokeProbs")
dir.result.meta.analysis.ad <- file.path(dir.result.meta.analysis, "AD_vs_CN/")
dir.RNA_vs_DNAm <- file.path(dir.base,"analysis_results/RNA_vs_DNAm/withSmoke/")
#for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)



#-------------------------------------------------------------------------------
# Data 
#-------------------------------------------------------------------------------
load(file.path(dir.base, "datasets/Aux/ADNI_matched_rna_dnam_residuals_bySex_with_Smoke_probes_all_sig_cpgs_updated.rda"))
load(file.path(dir.base, "datasets/Aux/ADNI_matched_rna_dnam_residuals_bySex_with_Smoke_probes_interaction_sig_cpgs.rda"))

all(rownames(metadata.dnam.female) == rownames(metadata.exp.female))
all(rownames(metadata.dnam.male) == rownames(metadata.exp.male))
all(colnames(residuals.matched.exp.female) == colnames(residuals.matched.met.female))
all(colnames(residuals.matched.exp.male) == colnames(residuals.matched.met.male))


#-----------------------------------------------------------------------------
# gene expression
#-----------------------------------------------------------------------------
load(file.path(dir.base, "datasets/Aux/ADNI_matched_rna_dnam_male_female_with_Smoke_probes.rda"))

# Male
Affy_Plate <- ADNI_Gene_Expression_Metadata[7,-c(1:3)] %>% as.numeric()
Affy_Plate <- Affy_Plate[
  match(
    gsub(" ","0",formatC(stringr::str_extract(pattern = "^[0-9]*",colnames(expression.matrix.male)))), 
    stringr::str_extract(pattern = "[0-9]*$",ADNI_Gene_Expression_Metadata[3,-c(1:3)] %>% as.character()))
]

RIN <- ADNI_Gene_Expression_Metadata[6,-c(1:3)] %>% as.numeric()
RIN <- RIN[
  match(
    gsub(" ","0",formatC(stringr::str_extract(pattern = "^[0-9]*",colnames(expression.matrix.male)))), 
    stringr::str_extract(pattern = "[0-9]*$",ADNI_Gene_Expression_Metadata[3,-c(1:3)] %>% as.character()))
]
metadata.exp <- colData(adni.se.male)[,c("age_at_visit","DX"),drop = FALSE]
metadata.exp$Affy_Plate <- Affy_Plate
metadata.exp$RIN <- RIN


aux <- expression.matrix.male
rownames(aux) <- MethReg:::map_ensg_to_symbol(rownames(aux))
library(xCell)
xcell <- xCell::xCellAnalysis(aux)
xcell <- xcell[c("B-cells", "NKT", "CD4+ T-cells", "CD8+ T-cells", "Monocytes", "Neutrophils"),]

metadata.exp <- cbind(metadata.exp,t(xcell))
colnames(metadata.exp)[5] <- "B_cells"
colnames(metadata.exp)[7] <- "CD4_T_cells"
colnames(metadata.exp)[8] <- "CD8_T_cells"

all(paste0(adni.se.male$RID,"_",adni.se.male$COLPROT,"_",adni.se.male$VISCODE) == colnames(expression.matrix.male))

doParallel::registerDoParallel(5)
res.male <-   plyr::adply(expression.matrix.male,.margins = 1,.fun = function(exp){
  dat <- cbind(exp, metadata.exp)
  dat$DX <- factor(dat$DX,levels = c("CN","Dementia"))
  res <- lm(log2 (exp) ~ DX + Affy_Plate + age_at_visit + B_cells +  NKT +  CD4_T_cells +  CD8_T_cells + Monocytes + Neutrophils,data = dat) %>%
    summary %>% coefficients() %>% as.data.frame()
  res[2,]
},.progress = "time", .parallel = T)
res.male$gene_symbol <- MethReg:::map_ensg_to_symbol(res.male[[1]])
res.male$fdr <- p.adjust(res.male$`Pr(>|t|)`,"fdr")

# FEMALE

Affy_Plate <- ADNI_Gene_Expression_Metadata[7,-c(1:3)] %>% as.numeric()
Affy_Plate <- Affy_Plate[
  match(
    gsub(" ","0",formatC(stringr::str_extract(pattern = "^[0-9]*",colnames(expression.matrix.female)))), 
    stringr::str_extract(pattern = "[0-9]*$",ADNI_Gene_Expression_Metadata[3,-c(1:3)] %>% as.character()))
]

RIN <- ADNI_Gene_Expression_Metadata[6,-c(1:3)] %>% as.numeric()
RIN <- RIN[
  match(
    gsub(" ","0",formatC(stringr::str_extract(pattern = "^[0-9]*",colnames(expression.matrix.female)))), 
    stringr::str_extract(pattern = "[0-9]*$",ADNI_Gene_Expression_Metadata[3,-c(1:3)] %>% as.character()))
]
metadata.exp <- colData(adni.se.female)[,c("age_at_visit","DX"),drop = FALSE]
metadata.exp$Affy_Plate <- Affy_Plate
metadata.exp$RIN <- RIN


aux <- expression.matrix.female
rownames(aux) <- MethReg:::map_ensg_to_symbol(rownames(aux))
library(xCell)
xcell <- xCell::xCellAnalysis(aux)
xcell <- xcell[c("B-cells", "NKT", "CD4+ T-cells", "CD8+ T-cells", "Monocytes", "Neutrophils"),]

metadata.exp <- cbind(metadata.exp,t(xcell))
colnames(metadata.exp)[5] <- "B_cells"
colnames(metadata.exp)[7] <- "CD4_T_cells"
colnames(metadata.exp)[8] <- "CD8_T_cells"

all(paste0(adni.se.female$RID,"_",adni.se.female$COLPROT,"_",adni.se.female$VISCODE) == colnames(expression.matrix.female))

doParallel::registerDoParallel(10)
res.female <-   plyr::adply(expression.matrix.female,.margins = 1,.fun = function(exp){
  dat <- cbind(exp, metadata.exp)
  dat$DX <- factor(dat$DX,levels = c("CN","Dementia"))
  res <- lm(log2 (exp) ~ DX + Affy_Plate + age_at_visit + B_cells +  NKT +  CD4_T_cells +  CD8_T_cells + Monocytes + Neutrophils,data = dat) %>%
    summary %>% coefficients() %>% as.data.frame()
  res[2,]
},.progress = "time",.parallel = T)
res.female$gene_symbol <- MethReg:::map_ensg_to_symbol(res.female[[1]])
res.female$fdr <- p.adjust(res.female$`Pr(>|t|)`,"fdr")


selected.genes.females <-readxl::read_excel("analysis_results/RNA_vs_DNAm/withSmoke/Supp Table 9 DNAm-RNA-Females-updated.xlsx")
selected.genes.females <- na.omit(unique(selected.genes.females$...3))[-c(1:2)]
selected.genes.males <- readxl::read_excel("analysis_results/RNA_vs_DNAm/withSmoke/Supp Table 9 DNAm-RNA-Males-updated.xlsx")
selected.genes.males <- na.omit(unique(selected.genes.males$...3))[-c(1:2)]

FEMALE_selected <- res.female[res.female$gene_symbol %in% selected.genes.females,]
FEMALE_selected$fdr <- p.adjust(FEMALE_selected$`Pr(>|t|)`,"fdr")
MALE_selected <- res.male[res.male$gene_symbol %in% selected.genes.males,]
MALE_selected$fdr <- p.adjust(MALE_selected$`Pr(>|t|)`,"fdr")



writexl::write_xlsx(
  list(
    "MALE_log2exp_vs_DX_plus_cofactor" = res.male,
    "MALE_selected" = MALE_selected,
    "FEMALE_log2exp_vs_DX_plus_cofactor" = res.female,
    "FEMALE_selected" = FEMALE_selected
  ),path = "analysis_results/RNA_vs_DNAm/withSmoke/Gene_expression_only.xlsx"
)

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
    
  } else if (analysis == "MALE"){
    idx.ad.vs.cn <- which(results$probeID %in% sig.cpgs.males)
    idx.cross <- which(results$probeID %in% sig.cross.cpgs.males)
    idx.combp <- which(results$probeID %in% sig.cpgs.in.dmr.males)
    
    n.cpgs.ad <- length(sig.cpgs.males)
    n.cpgs.cross <- length(sig.cross.cpgs.males)
    n.cpgs.compb <- length(sig.cpgs.in.dmr.males)
    
  } else {
    idx.interaction <- which(results$probeID %in% sig.cpgs.interaction)
    
    n.cpgs.interaction <- length(sig.cpgs.interaction)
  }
  
  if(analysis %in% c("FEMALE", "MALE")){
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
  } else {
    label.int <- paste0(analysis," - ",n.cpgs.interaction, " CpGs from interaction analysis")
    results <- results[idx.interaction,] 
    
    results$cpgs_from <- c(
      rep(label.int, length(idx.interaction))
    )
    
    results$RLM_met.residual_fdr <- c(
      p.adjust(results$RLM_met.residual_pvalue[which(results$cpgs_from %in% label.int)],method = "fdr")
    )
  }
  
  return(results)
}

#-----------------------------------------------------------------------------
# Select cpgs
#-----------------------------------------------------------------------------
# FEMALE CpGs with P<1E- 5 in AD vs. CN comparison
sig.cpgs.females <- readxl::read_xlsx("~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/DRAFT-FIGURES-TABLES-REVISION_7-1-2022/Supp Table 2 sig_female_male_cpgs_AD_vs_CN_updated_final_LWedit.xlsx",skip = 3,n_max = 24)
sig.cpgs.females <- sig.cpgs.females[-1,]
header <- colnames(sig.cpgs.females)
sig.cpgs.females <- sig.cpgs.females$cpg
length(sig.cpgs.females)


# MALE CpGs with P<1E- 5 in AD vs. CN comparison
sig.cpgs.males <- readxl::read_xlsx("~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/DRAFT-FIGURES-TABLES-REVISION_7-1-2022/Supp Table 2 sig_female_male_cpgs_AD_vs_CN_updated_final_LWedit.xlsx", skip = 29)
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

# Interaction results

sig.meta_df_AD_vs_CN.annotated.interaction <- readr::read_csv(
  file.path(
    dir.base,  "analysis_results/meta_analysis/Logistic_regression_model/withSmokePrbs/AD_vs_CN/meta_analysis_glm_fixed_effect_ADN_and_AIBL_AD_vs_CNI_interaction_annotated.csv"
  )
) %>% filter(pVal.final.bacon < 1E-05)
sig.cpgs.interaction <- sig.meta_df_AD_vs_CN.annotated.interaction$cpg


sig.cpgs.in.dmr.females <- intersect(sig.cpgs.in.dmr.females,meta_df_AD_vs_CN.annotated.female$cpg)
sig.cpgs.in.dmr.males <- intersect(sig.cpgs.in.dmr.males,meta_df_AD_vs_CN.annotated.male$cpg)

sig.all.cpgs.male <- unique(c(sig.cpgs.in.dmr.males,sig.cpgs.males,sig.cross.cpgs.males))
sig.all.cpgs.female <- unique(c(sig.cpgs.in.dmr.females,sig.cpgs.females,sig.cross.cpgs.females))

sig.single.cross.cpgs.male <- unique(c(sig.cpgs.males,sig.cross.cpgs.males))
sig.single.cross.cpgs.female <- unique(c(sig.cpgs.females,sig.cross.cpgs.females))

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



# 3 out of 4 cpgs in males are promoter
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

# 7 cpgs in males are distal
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

promoter.gene.dnam.pair.males <- promoter.gene.dnam.pair
distal.gene.dnam.pair.males <- distal.gene.dnam.pair
length(unique(promoter.gene.dnam.pair$target_symbol)) # 3 # 9
length(unique(distal.gene.dnam.pair$target_symbol)) # 56 # 153
length(union(unique(distal.gene.dnam.pair$target_symbol),unique(promoter.gene.dnam.pair$target_symbol)))
# 59 # 162

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


promoter.gene.dnam.pair[promoter.gene.dnam.pair$regionID %in% "chr6:31590640-31590641",]
# 13 out of in females are promoter
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
    auxfunction(
      row = row,
      exp = residuals.matched.exp.female,
      metadata.dnam = metadata.dnam.female,
      dnam = dnam
    )    
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

# 92 cpgs in are distal
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

promoter.gene.dnam.pair.females <- promoter.gene.dnam.pair
distal.gene.dnam.pair.females <- distal.gene.dnam.pair

length(unique(promoter.gene.dnam.pair$target_symbol)) # 8 # 30
length(unique(distal.gene.dnam.pair$target_symbol)) # 153 # 320
length(union(unique(distal.gene.dnam.pair$target_symbol),unique(promoter.gene.dnam.pair$target_symbol)))
# 161 # 350


#-------------------------------------------------------------------------------
# Interaction female
#-------------------------------------------------------------------------------
# Promoter analysis, triplets using remap
#-------------------------------------------------------------------------------
dnam <- MethReg:::map_probes_to_regions(
  dnam = residuals.matched.met.interaction.female,
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


# 13 out of in females are promoter
table((promoter.gene.dnam.pair$regionID %>% unique()) %in% make_names_from_granges(EPIC.hg19[sig.cpgs.interaction]))

nrow(promoter.gene.dnam.pair)
promoter.gene.dnam.pair <- promoter.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp.interaction.female))
nrow(promoter.gene.dnam.pair)
promoter.gene.dnam.pair$probeID <- names(EPIC.hg19)[match(promoter.gene.dnam.pair$regionID,make_names_from_granges(EPIC.hg19))]
promoter.gene.dnam.pair <- promoter.gene.dnam.pair %>% dplyr::filter(.data$regionID %in% rownames(dnam))
nrow(promoter.gene.dnam.pair)
results.promoter.analysis <- plyr::adply(
  .data = promoter.gene.dnam.pair,
  .margins = 1,
  .fun = function(row){
    auxfunction(
      row = row,
      exp = residuals.matched.exp.interaction.female,
      metadata.dnam = metadata.dnam.interaction.female,
      dnam = dnam
    )    
  }
)

results.promoter.analysis <- results.promoter.analysis %>% add_cpgs_from_and_do_fdr(analysis = "INTERACTION")

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

# 83 cpgs in are distal
table((distal.gene.dnam.pair$regionID %>% unique()) %in% make_names_from_granges(EPIC.hg19[sig.cpgs.interaction]))

nrow(distal.gene.dnam.pair)
distal.gene.dnam.pair <- distal.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp.interaction.female))
nrow(distal.gene.dnam.pair)
distal.gene.dnam.pair$probeID <- names(EPIC.hg19)[match(distal.gene.dnam.pair$regionID,make_names_from_granges(EPIC.hg19))]
distal.gene.dnam.pair <- distal.gene.dnam.pair %>% dplyr::filter(.data$regionID %in% rownames(dnam))
nrow(distal.gene.dnam.pair)

results.distal.analysis <- plyr::adply(
  .data = distal.gene.dnam.pair,
  .margins = 1,
  .fun = function(row){
    auxfunction(row = row, exp = residuals.matched.exp.interaction.female, metadata.dnam = metadata.dnam.interaction.female, dnam = dnam)    
  }
)


results.distal.analysis <- results.distal.analysis %>% add_cpgs_from_and_do_fdr("INTERACTION")

results.distal.analysis[results.distal.analysis$RLM_met.residual_fdr < 0.05,]

results.promoter.analysis.interaction.female <- results.promoter.analysis
results.distal.analysis.interaction.female <- results.distal.analysis

promoter.gene.dnam.pair.interaction.female <- promoter.gene.dnam.pair
distal.gene.dnam.pair.interaction.female <- distal.gene.dnam.pair

length(unique(promoter.gene.dnam.pair$target_symbol)) # 1
length(unique(distal.gene.dnam.pair$target_symbol)) # 22
length(union(unique(distal.gene.dnam.pair$target_symbol),unique(promoter.gene.dnam.pair$target_symbol)))
# 23
#-------------------------------------------------------------------------------
# Interaction male
#-------------------------------------------------------------------------------
# Promoter analysis, triplets using remap
#-------------------------------------------------------------------------------
dnam <- MethReg:::map_probes_to_regions(
  dnam = residuals.matched.met.interaction.male,
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


# 13 out of in females are promoter
table((promoter.gene.dnam.pair$regionID %>% unique()) %in% make_names_from_granges(EPIC.hg19[sig.cpgs.interaction]))

nrow(promoter.gene.dnam.pair)
promoter.gene.dnam.pair <- promoter.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp.interaction.male))
nrow(promoter.gene.dnam.pair)
promoter.gene.dnam.pair$probeID <- names(EPIC.hg19)[match(promoter.gene.dnam.pair$regionID,make_names_from_granges(EPIC.hg19))]
promoter.gene.dnam.pair <- promoter.gene.dnam.pair %>% dplyr::filter(.data$regionID %in% rownames(dnam))
nrow(promoter.gene.dnam.pair)
results.promoter.analysis <- plyr::adply(
  .data = promoter.gene.dnam.pair,
  .margins = 1,
  .fun = function(row){
    auxfunction(
      row = row,
      exp = residuals.matched.exp.interaction.male,
      metadata.dnam = metadata.dnam.interaction.male,
      dnam = dnam
    )    
  }
)

results.promoter.analysis <- results.promoter.analysis %>% add_cpgs_from_and_do_fdr(analysis = "INTERACTION")

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

# 83 cpgs in are distal
table((distal.gene.dnam.pair$regionID %>% unique()) %in% make_names_from_granges(EPIC.hg19[sig.cpgs.interaction]))

nrow(distal.gene.dnam.pair)
distal.gene.dnam.pair <- distal.gene.dnam.pair %>% dplyr::filter(.data$target %in% rownames(residuals.matched.exp.interaction.male))
nrow(distal.gene.dnam.pair)
distal.gene.dnam.pair$probeID <- names(EPIC.hg19)[match(distal.gene.dnam.pair$regionID,make_names_from_granges(EPIC.hg19))]
distal.gene.dnam.pair <- distal.gene.dnam.pair %>% dplyr::filter(.data$regionID %in% rownames(dnam))
nrow(distal.gene.dnam.pair)

results.distal.analysis <- plyr::adply(
  .data = distal.gene.dnam.pair,
  .margins = 1,
  .fun = function(row){
    auxfunction(row = row, exp = residuals.matched.exp.interaction.male, metadata.dnam = metadata.dnam.interaction.male, dnam = dnam)    
  }
)


results.distal.analysis <- results.distal.analysis %>% add_cpgs_from_and_do_fdr("INTERACTION")

results.distal.analysis[results.distal.analysis$RLM_met.residual_fdr < 0.05,]

results.promoter.analysis.interaction.male <- results.promoter.analysis
results.distal.analysis.interaction.male <- results.distal.analysis

promoter.gene.dnam.pair.interaction.male <- promoter.gene.dnam.pair
distal.gene.dnam.pair.interaction.male <- distal.gene.dnam.pair

length(unique(promoter.gene.dnam.pair$target_symbol)) # 1
length(unique(distal.gene.dnam.pair$target_symbol)) # 22
length(union(unique(distal.gene.dnam.pair$target_symbol),unique(promoter.gene.dnam.pair$target_symbol)))
# 23
#-------------------------------------------------------------------------------
# Save results
#-------------------------------------------------------------------------------

writexl::write_xlsx(
  list(
    "MALE_Promoter" = results.promoter.analysis.male,
    "MALE_Distal_10_up_10_down" = results.distal.analysis.male,
    "FEMALE_Promoter" = results.promoter.analysis.female,
    "FEMALE_Distal_10_up_10_down" = results.distal.analysis.female
    #"INTER_Promoter" = results.promoter.analysis.interaction,
    #"INTER_Distal_10_up_10_down" = results.distal.analysis.interaction
  ),
  path = file.path(dir.RNA_vs_DNAm,"FEMALE_MALE_Blood_ADNI_cpg_Target_vs_DNAm_updated.xlsx")
)

writexl::write_xlsx(
  list(
    "INTER_M_Promoter" = results.promoter.analysis.interaction.male,
    "INTER_M_Distal_10_up_10_down" = results.distal.analysis.interaction.male,
    "INTER_F_Promoter" = results.promoter.analysis.interaction.female,
    "INTER_F_Distal_10_up_10_down" = results.distal.analysis.interaction.female
    #"INTER_Promoter" = results.promoter.analysis.interaction,
    #"INTER_Distal_10_up_10_down" = results.distal.analysis.interaction
  ),
  path = file.path(dir.RNA_vs_DNAm,"INTERACTION_SEP_BY_FEMALE_MALE_Blood_ADNI_cpg_Target_vs_DNAm_updated.xlsx")
)

