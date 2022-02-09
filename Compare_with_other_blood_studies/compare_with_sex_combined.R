library(dplyr)
library(SummarizedExperiment)
devtools::load_all("~/Documents/packages/coMethDMR/")

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Combined brain sex analysis
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Supp Table https://www.dropbox.com/s/mf26sbblw14fyxm/ALL_SupplementaryTables_4-8-2021.xlsx?dl=0
# Put cpgs and dmrs in Supp Table S3-S6 together and then compare with each highlighted table above

all.brain.s3 <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD_metaAnalysis_bySex/DRAFT-REVISION_4-5-2021/ALL_SupplementaryTables_4-8-2021.xlsx",sheet = 3,skip = 3)
all.brain.s4 <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD_metaAnalysis_bySex/DRAFT-REVISION_4-5-2021/ALL_SupplementaryTables_4-8-2021.xlsx",sheet = 4,skip = 3)
all.brain.s5 <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD_metaAnalysis_bySex/DRAFT-REVISION_4-5-2021/ALL_SupplementaryTables_4-8-2021.xlsx",sheet = 5,skip = 3)
all.brain.s6 <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD_metaAnalysis_bySex/DRAFT-REVISION_4-5-2021/ALL_SupplementaryTables_4-8-2021.xlsx",sheet = 6,skip = 3)

all.brain.cpgs.dmrs <- c(all.brain.s3$cpg,all.brain.s4$cpg,all.brain.s5$DMR,all.brain.s6$DMR)
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Combined blood sex analysis
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Please put cpgs and dmrs in 1a-1c together and then compare with each highlighted table above
# (a) CpGs in Supp Tables S2, S4 from this file 
# https://www.dropbox.com/s/45xlp9jy6g55smn/_ALL_Supp_Tables-revision_1-17-2022.xlsx?dl=1
all.blood.sup_2 <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/DRAFT-REVISION_1-20-2022-submitted/_ALL_Supp_Tables-revision_1-17-2022.xlsx",sheet = 2,skip = 3)
all.blood.sup_4 <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/DRAFT-REVISION_1-20-2022-submitted/_ALL_Supp_Tables-revision_1-17-2022.xlsx",sheet = 4,skip = 3)

# (b)AD vs. CN comparison DMRs  
# https://www.dropbox.com/s/219bz8umm8c772x/_Main%20Table%202%20DMRs-Combp-AD_vs_CN_annotated.xlsx?dl=0
all.blood.ad.cn <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/DRAFT-TABLES_FIGURES_12-11-2021/_Main Table 2 DMRs-Combp-AD_vs_CN_annotated.xlsx",sheet = 1,skip = 1)

# (c) cross-tissue DMRs 
# https://www.dropbox.com/s/sihv4qqq5qt5g0l/_Main%20Table%203%20Top%2010%20prioritized-CpGs_and_DMRs-crossTissue_brain_blood-V2.xlsx?dl=0
all.blood.cross.analysis <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples/DRAFT-TABLES_FIGURES_12-11-2021/_Main Table 3 Top 10 prioritized-CpGs_and_DMRs-crossTissue_brain_blood-V2.xlsx",sheet = 1,skip = 3)

all.blood.cross.cpgs <- all.blood.cross.analysis[c(3:12),]
all.blood.cross.dmr <- all.blood.cross.analysis[c(15:24),]

all.blood.cpgs.dmrs <- na.omit(c(all.blood.sup_2$cpg,all.blood.sup_4$CpG,all.blood.ad.cn$DMR,all.blood.cross.cpgs$`CpG or DMR`,all.blood.cross.dmr$`CpG or DMR`))


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# (3) sex by aging study - McCartney et al. (2020) 
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
download.file("https://static-content.springer.com/esm/art%3A10.1186%2Fs13073-019-0693-z/MediaObjects/13073_2019_693_MOESM1_ESM.xls",destfile = "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/datasets/Aux/13073_2019_693_MOESM1_ESM.xls")
McCartney.S6 <- readxl::read_excel("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/datasets/Aux/13073_2019_693_MOESM1_ESM.xls",sheet = 6,skip = 1)$Probe
McCartney.S4 <- readxl::read_excel("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/datasets/Aux/13073_2019_693_MOESM1_ESM.xls",sheet = 4,skip = 1) %>% 
  dplyr::filter(Discovery.P.Value < 3.6E-8  & Replication.P.Value < 3.6E-8) %>% dplyr::pull(1)
McCartney.S5 <- readxl::read_excel("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/datasets/Aux/13073_2019_693_MOESM1_ESM.xls",sheet = 5,skip = 1) %>% 
  dplyr::filter(Discovery.P.Value < 3.6E-8  & Replication.P.Value < 3.6E-8)  %>% dplyr::pull(1)

McCartney <- c(McCartney.S6,McCartney.S4,McCartney.S5)
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# (4) sex-by-aging study: Yusipov et al. (2020) 
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
download.file("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7762479/bin/aging-12-202251-s005.xlsx",destfile = "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/datasets/Aux/aging-12-202251-s005.xlsx")
Yusipov <- readxl::read_excel("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/datasets/Aux/aging-12-202251-s005.xlsx")[[1]]

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Sex-specific analysis
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sex.sepecific.sup_tab2 <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/DRAFT-FIGURES-TABLES_1-11-2021/Supp Table 2 sig_female_male_cpgs_AD_vs_CN.xlsx",skip = 3)
sex.sepecific.sup_tab3 <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/DRAFT-FIGURES-TABLES_1-11-2021/Supp Table 3 Females-all-DMRs cnew.regions-p.bed_annotated.xlsx",skip = 1)
sex.sepecific.sup_tab4 <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/DRAFT-FIGURES-TABLES_1-11-2021/Supp Table 4 Males-all-DMRs cnew.regions-p.bed_annotated.xlsx",skip = 1)
sex.sepecific.tab3.female <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/DRAFT-FIGURES-TABLES_1-11-2021/_Main Table 3 cross-tissue-analysis.xlsx",skip = 3)
sex.sepecific.tab3.female <- sex.sepecific.tab3.female[1:12,]
sex.sepecific.tab3.male <- readxl::read_xlsx("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex/DRAFT-FIGURES-TABLES_1-11-2021/_Main Table 3 cross-tissue-analysis.xlsx",skip = 19)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Overlap
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# Sup tab 2
sex.sepecific.sup_tab2$CpG[sex.sepecific.sup_tab2$CpG %in% all.blood.cpgs.dmrs]
sex.sepecific.sup_tab2$CpG[sex.sepecific.sup_tab2$CpG %in% all.brain.cpgs.dmrs]
sex.sepecific.sup_tab2$CpG[sex.sepecific.sup_tab2$CpG %in% Yusipov]
sex.sepecific.sup_tab2$CpG[sex.sepecific.sup_tab2$CpG %in% McCartney]


# Sup tab 3
sex.sepecific.sup_tab3$DMR[sex.sepecific.sup_tab3$DMR %in% all.blood.cpgs.dmrs]
sex.sepecific.sup_tab3$DMR[sex.sepecific.sup_tab3$DMR %in% all.brain.cpgs.dmrs]

subsetByOverlaps(sex.sepecific.sup_tab3$DMR  %>% MethReg::make_granges_from_names(),
             grep("cg",all.blood.cpgs.dmrs,invert = T,value = T) %>% MethReg::make_granges_from_names()
)
subsetByOverlaps(sex.sepecific.sup_tab3$DMR  %>% MethReg::make_granges_from_names(),
             grep("cg",all.brain.cpgs.dmrs,invert = T,value = T) %>% MethReg::make_granges_from_names()
)


# Sup tab 4
sex.sepecific.sup_tab4$DMR[sex.sepecific.sup_tab4$DMR %in% all.blood.cpgs.dmrs]
sex.sepecific.sup_tab4$DMR[sex.sepecific.sup_tab4$DMR %in% all.brain.cpgs.dmrs]

subsetByOverlaps(sex.sepecific.sup_tab4$DMR  %>% MethReg::make_granges_from_names(),
                 grep("cg",all.blood.cpgs.dmrs,invert = T,value = T) %>% MethReg::make_granges_from_names()
)
subsetByOverlaps(sex.sepecific.sup_tab4$DMR  %>% MethReg::make_granges_from_names(),
                 grep("cg",all.brain.cpgs.dmrs,invert = T,value = T) %>% MethReg::make_granges_from_names()
)

# tab 3


sex.sepecific.tab3.male$CpG[sex.sepecific.tab3.male$CpG %in% all.blood.cpgs.dmrs]
sex.sepecific.tab3.male$CpG[sex.sepecific.tab3.male$CpG %in% all.brain.cpgs.dmrs]
sex.sepecific.tab3.male$CpG[sex.sepecific.tab3.male$CpG %in% Yusipov]
sex.sepecific.tab3.male$CpG[sex.sepecific.tab3.male$CpG %in% McCartney]


sex.sepecific.tab3.female$CpG[sex.sepecific.tab3.female$CpG %in% all.blood.cpgs.dmrs]
sex.sepecific.tab3.female$CpG[sex.sepecific.tab3.female$CpG %in% all.brain.cpgs.dmrs]
sex.sepecific.tab3.female$CpG[sex.sepecific.tab3.female$CpG %in% Yusipov]
sex.sepecific.tab3.female$CpG[sex.sepecific.tab3.female$CpG %in% McCartney]


