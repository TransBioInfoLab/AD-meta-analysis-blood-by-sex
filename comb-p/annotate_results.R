setwd("~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex")

# 1. Get comb-p results for female and males
combp.female <- readxl::read_xlsx("Michale/female_dist750/cnew.regions-p.bed.xlsx")
combp.female <- combp.female[combp.female$z_sidak_p < 0.05 & combp.female$n_probes >= 3,]
combp.male <- readxl::read_xlsx("Michale/male_dist750//cnew.regions-p.bed.xlsx")
combp.male <- combp.male[combp.male$z_sidak_p < 0.05 & combp.male$n_probes >= 3,]

dim(combp.male)
dim(combp.female)

combp.male$`#chrom`[combp.male$`#chrom` ==  0] <- "X"
plyr::count(combp.male$`#chrom`)
plyr::count(combp.female$`#chrom`)


# 2. Meta-analysis results for female and males

dir.base <- "~/TBL Dropbox/Tiago Silva/AD-meta-analysis-blood-samples-bySex"
dir.results <- file.path(dir.base,"analysis_results")
dir.result.meta.analysis <- file.path(dir.results, "meta_analysis/Logistic_regression_model/")
dir.result.meta.analysis.ad <- file.path(dir.result.meta.analysis, "AD_vs_CN/")

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



# 3.1 Data for annotation of enhancers
dir.base <- "~/TBL Dropbox/Tiago Silva//AD-meta-analysis-blood-samples/"
dir.data.aux <- file.path(dir.base,"datasets/Aux/") 


library(GenomicRanges)
library(dplyr)
data <- readr::read_tsv(
  file.path(dir.base,"datasets/nasser_2021/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz")
)
CellType.selected <- readxl::read_xlsx(
  file.path(dir.base,"code/annotations/Nassser study selected biosamples.xlsx"),col_names = FALSE
) %>% dplyr::pull(1)

data.filtered <- data %>% dplyr::filter(CellType %in% CellType.selected) %>% 
  dplyr::filter(!isSelfPromoter)  %>% 
  dplyr::filter(class != "promoter")

nasser.enhancer.gr <- data.filtered %>% makeGRangesFromDataFrame(
  start.field = "start",
  end.field = "end",
  seqnames.field = "chr",
  keep.extra.columns = TRUE
)


# 3.2 Data for annotation of enhancers
library(rGREAT)

devtools::load_all("~/Documents/packages/coMethDMR/")
dir.base <- "~/TBL Dropbox/Tiago Silva//AD-meta-analysis-blood-samples/"
dir.data.aux <- file.path(dir.base,"datasets/Aux/") 
load(file.path(dir.data.aux,"E073_15_coreMarks_segments.rda"))


# 4 Annotation of enhancers

add_annotation <- function(result){
  
  message("Annotating E073_15_coreMarks_segments")
  result$seqnames <- paste0("chr",result$`#chrom`) 
  result$region <- paste0(result$seqnames,":",result$start,"-", result$end)      
  result$start <- as.numeric(result$start)
  result$end <- as.numeric(result$end)
  result.gr <- result %>% makeGRangesFromDataFrame(
    start.field = "start",
    end.field = "end",
    seqnames.field = "seqnames"
  )
  hits <- findOverlaps(result.gr, ChmmModels.gr) %>% as.data.frame()
  hits$state <- ChmmModels.gr$state[hits$subjectHits]
  hits$region <- result$region[hits$queryHits]
  result$E073_15_coreMarks_segments_state <- hits$state[match(result$region,hits$region)]
  
  message("Annotating GREAT")
  
  
  job <- submitGreatJob(result.gr, species = "hg19")
  regionsToGenes <- data.frame(plotRegionGeneAssociationGraphs(job))
  
  regionsToGenes$GREAT_annotation <- ifelse(
    regionsToGenes$distTSS > 0,
    paste0(regionsToGenes$gene, " (+", regionsToGenes$distTSS, ")"),
    paste0(regionsToGenes$gene, " (", regionsToGenes$distTSS, ")"))
  regionsToGenes <- regionsToGenes[
    ,c("seqnames", "start", "end", "GREAT_annotation")
  ]
  great <- regionsToGenes %>%
    group_by(seqnames, start, end) %>%
    mutate(GREAT_annotation = paste0(GREAT_annotation,collapse = ";")) %>% unique()
  
  result <- dplyr::left_join(result, great,by = c("seqnames","start","end"))
  
  message("Annotating Island")
  
  result$chrom <- result$seqnames
  result <- coMethDMR::AnnotateResults(result, arrayType = "EPIC", nCores_int = 10)
  result$chrom <-  result$chr <- NULL
  
  message("Annotating enhancer")
  
  result.gr <- result %>% makeGRangesFromDataFrame(
    start.field = "start",
    end.field = "end",
    seqnames.field = "seqnames"
  )
  hits <- findOverlaps(result.gr, nasser.enhancer.gr) %>% as.data.frame()
  result$nasser_is_enahncer <- FALSE
  result$nasser_is_enahncer[unique(hits$queryHits)] <- TRUE
  result$nasser_is_enahncer_cell_types <- NA
  result$nasser_is_enahncer_cell_types[unique(hits$queryHits)]  <- sapply(unique(hits$queryHits), function(x){nasser.enhancer.gr$CellType[hits$subjectHits[hits$queryHits %in% x]] %>% unique %>% paste(collapse = ",")})
  
  
  if("fdr.bacon" %in% colnames(result)) result <- result[order(result$fdr.bacon),]
  if("AdjP(BH)" %in% colnames(result)) result <- result[order(result$`AdjP(BH)`),]
  
  result$cpgs_in_region <- sapply(
    result$region, 
    function(x){
      coMethDMR:::GetCpGsInRegion(x,genome = "hg19",arrayType = "EPIC") %>% intersect(meta_df_AD_vs_CN.annotated.male$cpg) %>% paste(collapse = ",")
    })
  return(result)
  
}




combp.male.annotated <- combp.male %>% add_annotation

writexl::write_xlsx(
  x = combp.male.annotated,
  path = file.path("Michale/male_dist750/cnew.regions-p.bed_annotated.xlsx")
)

combp.male.annotated.top_10 <- combp.male.annotated %>% dplyr::top_n(10,wt = -z_sidak_p)


tab <- plyr::llply(combp.male.annotated.top_10$cpgs_in_region,.fun = function(x){
  cpgs <- stringr::str_split(x,",") %>% unlist 
  meta_df_AD_vs_CN.annotated.male %>% dplyr::filter(cpg %in% cpgs)
})
names(tab) <- gsub(":|-","_",combp.male.annotated.top_10$region)


writexl::write_xlsx(
  x = append(list("combp.male.annotated.top_10" = combp.male.annotated.top_10),tab),
  path = file.path("Michale/male_dist750/cnew.regions-p.bed_annotated_top_10.xlsx")
)


dim(combp.male)
dim(combp.male.annotated)
dim(combp.male.annotated.top_10)



combp.female.annotated <- combp.female %>% add_annotation
writexl::write_xlsx(
  x = combp.female.annotated,
  path = file.path("Michale/female_dist750/cnew.regions-p.bed_annotated.xlsx")
)

combp.female.annotated.top_10 <- combp.female.annotated  %>% dplyr::top_n(10,wt = -z_sidak_p)

tab.female <- plyr::llply(combp.female.annotated.top_10$cpgs_in_region,.fun = function(x){
  cpgs <- stringr::str_split(x,",") %>% unlist 
  meta_df_AD_vs_CN.annotated.female %>% dplyr::filter(cpg %in% cpgs)
})
names(tab.female) <- gsub(":|-","_",combp.female.annotated.top_10$region)


writexl::write_xlsx(
  x = append(list("combp.female.annotated.top_10" = combp.female.annotated.top_10),tab.female),
  path = file.path("Michale/female_dist750/cnew.regions-p.bed_annotated_top_10.xlsx")
)

dim(combp.female)
dim(combp.female.annotated)
dim(combp.female.annotated.top_10)


# (3) make a venn diagram of overlapping DMRs from males and females
combp.female.annotated <- readxl::read_xlsx("Michale/female_dist750/cnew.regions-p.bed_annotated.xlsx")
combp.male.annotated <- readxl::read_xlsx("Michale/male_dist750/cnew.regions-p.bed_annotated.xlsx")

gplots::venn(list("Comp-b Males" = combp.male.annotated$region,"Comp-b females" = combp.female.annotated$region))

findOverlaps(
  combp.male.annotated$region %>% MethReg::make_granges_from_names(),
  combp.female.annotated$region %>% MethReg::make_granges_from_names()
)

findOverlaps(
  combp.male.annotated$region %>% MethReg::make_granges_from_names(),
  combp.female.annotated$region %>% MethReg::make_granges_from_names()
)

# (4) venn diagram for overlapping DMRs with single CpGs in file
sig.cpgs.females <- readxl::read_xlsx("DRAFT-FIGURES-TABLES_8-2-2021/Supp Table 2 sig_female_male_cpgs_AD_vs_CN.xlsx",skip = 3,n_max = 25)
sig.cpgs.females <- sig.cpgs.females[-1,]
sig.cpgs.males <- readxl::read_xlsx("DRAFT-FIGURES-TABLES_8-2-2021/Supp Table 2 sig_female_male_cpgs_AD_vs_CN.xlsx", skip = 30)
colnames(sig.cpgs.males) <- colnames(sig.cpgs.females)

findOverlaps(
  paste0(sig.cpgs.females$chr,":",sig.cpgs.females$pos,"-",sig.cpgs.females$pos) %>% MethReg::make_granges_from_names(),
  combp.female.annotated$region %>% MethReg::make_granges_from_names()
)

findOverlaps(
  paste0(sig.cpgs.males$chr,":",sig.cpgs.males$pos,"-",sig.cpgs.males$pos) %>% MethReg::make_granges_from_names(),
  combp.male.annotated$region %>% MethReg::make_granges_from_names()
)

ggpubr::ggarrange(
  plotlist = list(
    ggvenn::ggvenn(
      list(
        "Comp-b analysis \n DMR sig. in males" = combp.male.annotated$region,
        "Comp-b analysis \n DMR sig. in females" = combp.female.annotated$region),
      fill_color = c("#0073C2FF","#CD534CFF"),
      set_name_size = 3
    ),
    ggpubr::ggarrange(
      plotlist = list(
        
        ggvenn::ggvenn(
          list(
            "Males \n CpGs within Sig. DMR in\n comb-p analysis " = intersect(meta_df_AD_vs_CN.annotated.male$cpg,combp.male.annotated$cpgs_in_region %>% stringr::str_split(",") %>% unlist %>% unique),
            "Males \n Sig CpGs \nin single cpg analysis" = sig.cpgs.males$CpG
          ),
          fill_color = c("#EFC000FF","#868686FF"),
          set_name_size = 3
        ),
        
        ggvenn::ggvenn(
          list(
            "Females \nCpGs within Sig. DMR in\n comb-p analysis " = intersect(meta_df_AD_vs_CN.annotated.female$cpg,combp.female.annotated$cpgs_in_region %>% stringr::str_split(",") %>% unlist %>% unique),
            "Females \nSig CpGs in\n single cpg analysis" = sig.cpgs.females$CpG
          ),
          fill_color = c("#EFC000FF","#868686FF"),
          
          set_name_size = 3
        )
      ),ncol = 2
    )
  ), nrow = 2
)


gplots::venn(list(sig.cpgs.females$CpG,sig.cpgs.males$CpG))

