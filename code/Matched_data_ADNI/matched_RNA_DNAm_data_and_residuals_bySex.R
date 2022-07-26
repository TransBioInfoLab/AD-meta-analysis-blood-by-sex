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
# Date: 12 July 2021
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
library(dplyr)
library(matrixStats)
library(SummarizedExperiment)
library(MethReg)
library(xCell)
library(readr)
library(dorothea)


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
   

#-----------------------------------------------------------------------------
# Directories
#-----------------------------------------------------------------------------
cohort <- "ADNI"
dir.base <- "~/TBL Dropbox/Wei Zhang"
dir.base.data <- file.path(dir.base,"/AD-meta-analysis-blood-samples/")
dir.data <- file.path(dir.base.data,"datasets/",cohort,"/data/DNA_methylation") 
dir.data.clinical <- file.path(dir.base.data,"datasets/",cohort,"/data/Clinical") 
dir.data.raw <- file.path(dir.base.data,"datasets/",cohort,"/data/DNA_methylation/idat") 
dir.data.processed <- file.path(dir.base.data,"datasets/",cohort,"/data/DNA_methylation/processed/withSmoke") 
dir.results <- file.path(dir.base,"..//AD-meta-analysis-blood-samples-bySex/analysis_results/")
dir.result.meta.analysis <- file.path(dir.results, "meta_analysis/Logistic_regression_model/withSmokePrbs/")
dir.result.meta.analysis.ad <- file.path(dir.result.meta.analysis, "AD_vs_CN/")

#-----------------------------------------------------------------------------
# GET ADNI data
#-----------------------------------------------------------------------------
# DNA methylation
#-----------------------------------------------------------------------------
adni.se <- readRDS(file.path(dir.data.processed,"adni_se_min_age_at_visit_65_with_XY_AD_CN.RDS"))
 #adni.se <- adni.se[,adni.se$DX %in% c("CN","Dementia")]

# Gene expression: Affymetrix Human Genome U 219 array 
ADNI_Gene_Expression_Profile <- read_csv(file.path(dir.base.data,"datasets/ADNI/data/gene_expression/ADNI_Gene_Expression_Profile/ADNI_Gene_Expression_Profile.csv"),skip = 8)
ADNI_Gene_Expression_Profile$...748 <- NULL

# We have our gene expression as below:
#   ProbeSet      LocusLink Symbol                  X4    X5    X6    X7    X8    X9   X10
# 1 11715104_s_at LOC92736  OTOP2                 2.15  2.16  2.52  2.28  2.25  2.24  1.99
# 2 11715105_at   LOC284099 C17ORF78              2.27  2.13  1.96  2.35  2.15  2.06  2.32
# 3 11715106_x_at NA        CTAGE6 || CTAGE15     2.43  2.27  2.33  2.26  2.33  2.45  2.17
# The row (3) must be break in to two
#   ProbeSet      LocusLink Symbol                  X4    X5    X6    X7    X8    X9   X10
# 1 11715104_s_at LOC92736  OTOP2                 2.15  2.16  2.52  2.28  2.25  2.24  1.99
# 2 11715105_at   LOC284099 C17ORF78              2.27  2.13  1.96  2.35  2.15  2.06  2.32
# 3 11715106_x_at NA        CTAGE6                2.43  2.27  2.33  2.26  2.33  2.45  2.17
# 4 11715106_x_at NA        CTAGE15               2.43  2.27  2.33  2.26  2.33  2.45  2.17
ADNI_Gene_Expression_Profile <- ADNI_Gene_Expression_Profile %>% tidyr::separate_rows("Symbol")

probe.to.ensg <- readr::read_csv(file.path("~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/datasets/Aux/HG-U219.na36.annot.csv"), skip = 25)  %>% 
  tidyr::separate_rows("Ensembl")
#library(hgu219.db)
#x <- hgu219ENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID, this is required by MethReg
#probe.to.ensg <- as.data.frame(as.list(x) %>% unlist) %>% na.omit
ADNI_Gene_Expression_Profile$ENGS <- probe.to.ensg$Ensembl[match(ADNI_Gene_Expression_Profile$ProbeSet,probe.to.ensg$`Probe Set ID`)]
ADNI_Gene_Expression_Profile <- ADNI_Gene_Expression_Profile[!is.na(ADNI_Gene_Expression_Profile$ENGS),]

nrow(ADNI_Gene_Expression_Profile) # 52955 

# dropping genes in bottom 10 percentile for over 80% of the samples
genes.low.expressed.samples.count <- plyr::aaply(
  ADNI_Gene_Expression_Profile[,grep("\\.\\.",colnames(ADNI_Gene_Expression_Profile))] ,
  .margins = 2, # for each sample set get all genes
  .fun = function(sample.genes){
    # for each sample, mark the bottom 10% less expressed genes as 1
    sample.genes[[1]]  <= quantile(sample.genes[[1]] , probs = c(.10), type = 3)
  }) %>% colSums()
genes.idx <- which(genes.low.expressed.samples.count > (length(grep("\\.\\.",colnames(ADNI_Gene_Expression_Profile))) * 0.8))
ADNI_Gene_Expression_Profile <- ADNI_Gene_Expression_Profile[-c(genes.idx),]

nrow(ADNI_Gene_Expression_Profile) #  50303, previous 41681

# Since we have multiple probes mapping to the same gene in the array 
# we will take the median of the same genes
# as suggest in https://www.nature.com/articles/s41598-020-60595-1
expression.matrix <- plyr::aaply(unique(ADNI_Gene_Expression_Profile$ENGS),.margins = 1,.fun = function(gene){
  dat <- ADNI_Gene_Expression_Profile[ADNI_Gene_Expression_Profile$ENGS == gene,grep("\\.\\.",colnames(ADNI_Gene_Expression_Profile))]
  colMedians(dat %>% as.matrix)
},.progress = "time") 
rownames(expression.matrix) <- unique(ADNI_Gene_Expression_Profile$ENGS)

ADNI_Gene_Expression_Metadata <- read_csv(file.path(dir.base.data,"datasets/ADNI/data/gene_expression/ADNI_Gene_Expression_Profile/ADNI_Gene_Expression_Profile.csv"),skip = 0,col_names = FALSE,n_max = 7)
ADNI_Gene_Expression_Metadata$X748 <- NULL

gene.exp.IDs <- apply(
  ADNI_Gene_Expression_Metadata[,4:ncol(ADNI_Gene_Expression_Metadata)],MARGIN = 2, 
  FUN = function(col) {
    paste0(
      stringr::str_extract(pattern = "[0-9]*$",string = col[3]),
      "_", col[1],"_",
      col[2]
    )
  }
)

DXSUM <- readr::read_csv(file.path(dir.base.data,"datasets/ADNI/data/Clinical/DXSUM_PDXCONV_ADNIALL_downloaded_2-19-2021.csv"))
DXSUM <- DXSUM[match(gene.exp.IDs,paste0(stringr::str_extract(pattern = "[0-9]*$",DXSUM$PTID),"_",DXSUM$Phase,"_",DXSUM$VISCODE)),]
gene.exp.IDs <- paste0(DXSUM$RID,"_",DXSUM$Phase,"_",DXSUM$VISCODE2)
colnames(expression.matrix) <- gene.exp.IDs
expression.matrix <- expression.matrix[,colnames(expression.matrix) != "NA_NA_NA"]

#-----------------------------------------------------------------------------
# get matched DNAm and Gene expression
#-----------------------------------------------------------------------------
dnam.IDs <- paste0(adni.se$RID,"_",adni.se$COLPROT,"_",adni.se$VISCODE)
table(gene.exp.IDs %in% dnam.IDs) 
# FALSE  TRUE 
# 479    265 
common.ids <- base::intersect(dnam.IDs %>% as.character,gene.exp.IDs %>% as.character)
adni.se <- adni.se[,paste0(adni.se$RID,"_",adni.se$COLPROT,"_",adni.se$VISCODE) %in% common.ids]
adni.se <- adni.se[,match(common.ids,paste0(adni.se$RID,"_",adni.se$COLPROT,"_",adni.se$VISCODE))]
expression.matrix <- expression.matrix[,colnames(expression.matrix) %in% common.ids]
expression.matrix <- expression.matrix[,match(common.ids,colnames(expression.matrix))]

table(colnames(expression.matrix) == paste0(adni.se$RID,"_",adni.se$COLPROT,"_",adni.se$VISCODE))


idx.male <- which(adni.se$PTGENDER == "Male")
idx.female <- which(adni.se$PTGENDER == "Female")
adni.se.female <- adni.se[,idx.female]
expression.matrix.female <- expression.matrix[,idx.female]
adni.se.male <- adni.se[,idx.male]
expression.matrix.male <- expression.matrix[,idx.male]

save(
  adni.se,
  adni.se.female,
  adni.se.male,
  expression.matrix,
  expression.matrix.male,
  expression.matrix.female,
  ADNI_Gene_Expression_Metadata,
  file = "~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/datasets/Aux/ADNI_matched_rna_dnam_male_female_with_Smoke_probes.rda"
)

#-----------------------------------------------------------------------------
# Residuals 
#-----------------------------------------------------------------------------
# MALE 
#-----------------------------------------------------------------------------
library(readr)
se.selected.cpgs <- adni.se.male[rownames(adni.se.male) %in% sig.all.cpgs.male,]
colnames(se.selected.cpgs) <- paste0(se.selected.cpgs$RID,"_",se.selected.cpgs$COLPROT,"_",se.selected.cpgs$VISCODE)
se.selected.cpgs$granulocytes <- se.selected.cpgs$Eosino + se.selected.cpgs$Neutro
metadata.dnam <- colData(se.selected.cpgs)[,c("RID","DX","Mono","granulocytes","CD4T","NK","B","age_at_visit","PTGENDER","PlateNumber")]

residuals.matched.met <- MethReg::get_residuals(
  data = log2(assay(se.selected.cpgs) / (1 - assay(se.selected.cpgs))), # m-values
  metadata.samples = metadata.dnam[,c("granulocytes","Mono","CD4T","NK","B","age_at_visit","PlateNumber")],
  cores = 4
)

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

metadata.exp <- colData(se.selected.cpgs)[,c("age_at_visit"),drop = FALSE]
metadata.exp$Affy_Plate <- Affy_Plate
metadata.exp$RIN <- RIN


aux <- expression.matrix.male
rownames(aux) <- MethReg:::map_ensg_to_symbol(rownames(aux))
xcell <- xCellAnalysis(aux)
xcell <- xcell[c("B-cells", "NKT", "CD4+ T-cells", "CD8+ T-cells", "Monocytes", "Neutrophils"),]

metadata.exp <- cbind(metadata.exp,t(xcell))
colnames(metadata.exp)[4] <- "B_cells"
colnames(metadata.exp)[6] <- "CD4_T_cells"
colnames(metadata.exp)[7] <- "CD8_T_cells"

residuals.matched.exp <- MethReg::get_residuals(
  data = log2(expression.matrix.male),
  metadata.samples = metadata.exp,
  cores = 5
)

all(rownames(metadata.dnam) == rownames(metadata.exp))
all(colnames(residuals.matched.exp) == colnames(residuals.matched.met))

#-----------------------------------------------------------------------------
# get TF activity
#-----------------------------------------------------------------------------
regulons.dorothea <- dorothea::dorothea_hs
rnaseq.tf.es <- MethReg::get_tf_ES(
  exp = residuals.matched.exp,
  regulons = regulons.dorothea
)

metadata.dnam.male <- metadata.dnam
residuals.matched.exp.male <- residuals.matched.exp
metadata.exp.male <- metadata.exp
rnaseq.tf.es.male <- rnaseq.tf.es
residuals.matched.met.male <- residuals.matched.met
#-----------------------------------------------------------------------------
# Residuals 
#-----------------------------------------------------------------------------
# FEMALE 
#-----------------------------------------------------------------------------
se.selected.cpgs <- adni.se.female[rownames(adni.se.female) %in% sig.all.cpgs.female,]
colnames(se.selected.cpgs) <- paste0(se.selected.cpgs$RID,"_",se.selected.cpgs$COLPROT,"_",se.selected.cpgs$VISCODE)
se.selected.cpgs$granulocytes <- se.selected.cpgs$Eosino + se.selected.cpgs$Neutro
metadata.dnam <- colData(se.selected.cpgs)[,c("RID","DX","Mono","granulocytes","CD4T","NK","B","age_at_visit","PTGENDER","PlateNumber")]

residuals.matched.met <- MethReg::get_residuals(
  data = log2(assay(se.selected.cpgs) / (1 - assay(se.selected.cpgs))), # m-values
  metadata.samples = metadata.dnam[,c("granulocytes","Mono","CD4T","NK","B","age_at_visit","PlateNumber")],
  cores = 4
)

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

metadata.exp <- colData(se.selected.cpgs)[,c("age_at_visit"),drop = FALSE]
metadata.exp$Affy_Plate <- Affy_Plate
metadata.exp$RIN <- RIN


aux <- expression.matrix.female
rownames(aux) <- MethReg:::map_ensg_to_symbol(rownames(aux))
xcell <- xCellAnalysis(aux)
xcell <- xcell[c("B-cells", "NKT", "CD4+ T-cells", "CD8+ T-cells", "Monocytes", "Neutrophils"),]

metadata.exp <- cbind(metadata.exp,t(xcell))
colnames(metadata.exp)[4] <- "B_cells"
colnames(metadata.exp)[6] <- "CD4_T_cells"
colnames(metadata.exp)[7] <- "CD8_T_cells"

residuals.matched.exp <- MethReg::get_residuals(
  data = log2(expression.matrix.female),
  metadata.samples = metadata.exp,
  cores = 5
)

all(rownames(metadata.dnam) == rownames(metadata.exp))
all(colnames(residuals.matched.exp) == colnames(residuals.matched.met))

#-----------------------------------------------------------------------------
# get TF activity
#-----------------------------------------------------------------------------
regulons.dorothea <- dorothea::dorothea_hs
rnaseq.tf.es <- MethReg::get_tf_ES(
  exp = residuals.matched.exp,
  regulons = regulons.dorothea
)

metadata.dnam.female <- metadata.dnam
residuals.matched.exp.female <- residuals.matched.exp
metadata.exp.female <- metadata.exp
rnaseq.tf.es.female <- rnaseq.tf.es
residuals.matched.met.female <- residuals.matched.met
#-----------------------------------------------------------------------------
# Residuals
#-----------------------------------------------------------------------------
# Interaction analysis female
#-----------------------------------------------------------------------------
se.selected.cpgs <- adni.se.female[rownames(adni.se.female) %in% sig.cpgs.interaction,]
colnames(se.selected.cpgs) <- paste0(se.selected.cpgs$RID,"_",se.selected.cpgs$COLPROT,"_",se.selected.cpgs$VISCODE)
se.selected.cpgs$granulocytes <- se.selected.cpgs$Eosino + se.selected.cpgs$Neutro
metadata.dnam <- colData(se.selected.cpgs)[,c("RID","DX","Mono","granulocytes","CD4T","NK","B","age_at_visit","PTGENDER","PlateNumber")]

residuals.matched.met <- MethReg::get_residuals(
  data = log2(assay(se.selected.cpgs) / (1 - assay(se.selected.cpgs))), # m-values
  metadata.samples = metadata.dnam[,c("granulocytes","Mono","CD4T","NK","B","age_at_visit","PlateNumber")],
  cores = 4
)

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

metadata.exp <- colData(se.selected.cpgs)[,c("age_at_visit"),drop = FALSE]
metadata.exp$Affy_Plate <- Affy_Plate
metadata.exp$RIN <- RIN


aux <- expression.matrix.female
rownames(aux) <- MethReg:::map_ensg_to_symbol(rownames(aux))
xcell <- xCellAnalysis(aux)
xcell <- xcell[c("B-cells", "NKT", "CD4+ T-cells", "CD8+ T-cells", "Monocytes", "Neutrophils"),]

metadata.exp <- cbind(metadata.exp,t(xcell))
colnames(metadata.exp)[4] <- "B_cells"
colnames(metadata.exp)[6] <- "CD4_T_cells"
colnames(metadata.exp)[7] <- "CD8_T_cells"

residuals.matched.exp <- MethReg::get_residuals(
  data = log2(expression.matrix.female),
  metadata.samples = metadata.exp,
  cores = 5
)

all(rownames(metadata.dnam) == rownames(metadata.exp))
all(colnames(residuals.matched.exp) == colnames(residuals.matched.met))

#-----------------------------------------------------------------------------
# get TF activity
#-----------------------------------------------------------------------------
regulons.dorothea <- dorothea::dorothea_hs
rnaseq.tf.es <- MethReg::get_tf_ES(
  exp = residuals.matched.exp,
  regulons = regulons.dorothea
)

metadata.dnam.interaction.female <- metadata.dnam
residuals.matched.exp.interaction.female <- residuals.matched.exp
metadata.exp.interaction.female <- metadata.exp
rnaseq.tf.es.interaction.female <- rnaseq.tf.es
residuals.matched.met.interaction.female <- residuals.matched.met

#-----------------------------------------------------------------------------
# Residuals
#-----------------------------------------------------------------------------
# Interaction analysis male
#-----------------------------------------------------------------------------
se.selected.cpgs <- adni.se.male[rownames(adni.se.male) %in% sig.cpgs.interaction,]
colnames(se.selected.cpgs) <- paste0(se.selected.cpgs$RID,"_",se.selected.cpgs$COLPROT,"_",se.selected.cpgs$VISCODE)
se.selected.cpgs$granulocytes <- se.selected.cpgs$Eosino + se.selected.cpgs$Neutro
metadata.dnam <- colData(se.selected.cpgs)[,c("RID","DX","Mono","granulocytes","CD4T","NK","B","age_at_visit","PTGENDER","PlateNumber")]

residuals.matched.met <- MethReg::get_residuals(
  data = log2(assay(se.selected.cpgs) / (1 - assay(se.selected.cpgs))), # m-values
  metadata.samples = metadata.dnam[,c("granulocytes","Mono","CD4T","NK","B","age_at_visit","PlateNumber")],
  cores = 4
)

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

metadata.exp <- colData(se.selected.cpgs)[,c("age_at_visit"),drop = FALSE]
metadata.exp$Affy_Plate <- Affy_Plate
metadata.exp$RIN <- RIN


aux <- expression.matrix.male
rownames(aux) <- MethReg:::map_ensg_to_symbol(rownames(aux))
xcell <- xCellAnalysis(aux)
xcell <- xcell[c("B-cells", "NKT", "CD4+ T-cells", "CD8+ T-cells", "Monocytes", "Neutrophils"),]

metadata.exp <- cbind(metadata.exp,t(xcell))
colnames(metadata.exp)[4] <- "B_cells"
colnames(metadata.exp)[6] <- "CD4_T_cells"
colnames(metadata.exp)[7] <- "CD8_T_cells"

residuals.matched.exp <- MethReg::get_residuals(
  data = log2(expression.matrix.male),
  metadata.samples = metadata.exp,
  cores = 5
)

all(rownames(metadata.dnam) == rownames(metadata.exp))
all(colnames(residuals.matched.exp) == colnames(residuals.matched.met))

#-----------------------------------------------------------------------------
# get TF activity
#-----------------------------------------------------------------------------
regulons.dorothea <- dorothea::dorothea_hs
rnaseq.tf.es <- MethReg::get_tf_ES(
  exp = residuals.matched.exp,
  regulons = regulons.dorothea
)

metadata.dnam.interaction.male <- metadata.dnam
residuals.matched.exp.interaction.male <- residuals.matched.exp
metadata.exp.interaction.male <- metadata.exp
rnaseq.tf.es.interaction.male <- rnaseq.tf.es
residuals.matched.met.interaction.male <- residuals.matched.met


#-----------------------------------------------------------------------------
# Save data
#-----------------------------------------------------------------------------

save(
  metadata.dnam.female,
  metadata.dnam.male,
  residuals.matched.met.female,
  residuals.matched.met.male,
  metadata.exp.female,
  metadata.exp.male,
  rnaseq.tf.es.female,
  rnaseq.tf.es.male,
  residuals.matched.exp.female,
  residuals.matched.exp.male,
  file = "~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/datasets/Aux/ADNI_matched_rna_dnam_residuals_bySex_with_Smoke_probes_all_sig_cpgs_updated.rda"
)

save(
  metadata.dnam.interaction.female,
  metadata.dnam.interaction.male,
  residuals.matched.exp.interaction.female,
  residuals.matched.exp.interaction.male,
  metadata.exp.interaction.female,
  metadata.exp.interaction.male,
  rnaseq.tf.es.interaction.female,
  rnaseq.tf.es.interaction.male,
  residuals.matched.met.interaction.female,
  residuals.matched.met.interaction.male,
  file = "~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/datasets/Aux/ADNI_matched_rna_dnam_residuals_bySex_with_Smoke_probes_interaction_sig_cpgs.rda"
)
