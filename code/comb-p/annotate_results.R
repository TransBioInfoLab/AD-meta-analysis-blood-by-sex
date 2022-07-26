setwd("~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex")

# 1. Get comb-p results for female and males
combp.female <- readxl::read_xlsx("Michale/female_750/cnew.regions-p.bed.xlsx")
combp.female <- combp.female[combp.female$z_sidak_p < 0.05 & combp.female$n_probes >= 3,]
combp.male <- readxl::read_xlsx("Michale/male_750//cnew.regions-p.bed.xlsx")
combp.male <- combp.male[combp.male$z_sidak_p < 0.05 & combp.male$n_probes >= 3,]

dim(combp.male)
dim(combp.female)

combp.male$`#chrom`[combp.male$`#chrom` ==  0] <- "X"
plyr::count(combp.male$`#chrom`)
plyr::count(combp.female$`#chrom`)


# 2. Meta-analysis results for female and males

dir.base <- "~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex"
dir.results <- file.path(dir.base,"analysis_results")
dir.result.meta.analysis <- file.path(dir.results, "meta_analysis/Logistic_regression_model/withSmokePrbs")
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
dir.base <- "~/TBL Dropbox/Wei Zhang//AD-meta-analysis-blood-samples/"
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
library(coMethDMR)

#devtools::load_all("~/Documents/packages/coMethDMR/")
dir.base <- "~/TBL Dropbox/Wei Zhang//AD-meta-analysis-blood-samples/"
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

combp.male.annotated <- combp.male.annotated %>% mutate(
  n_probes = str_count(cpgs_in_region,"cg") 
)

writexl::write_xlsx(
  x = combp.male.annotated,
  path = file.path("Michale/male_750/cnew.regions-p.bed_annotated.xlsx")
)

combp.male.annotated.top_10 <- combp.male.annotated %>% dplyr::top_n(10,wt = -z_sidak_p)


tab <- plyr::llply(combp.male.annotated.top_10$cpgs_in_region,.fun = function(x){
  cpgs <- str_split(x,",") %>% unlist 
  meta_df_AD_vs_CN.annotated.male %>% dplyr::filter(cpg %in% cpgs)
})
names(tab) <- gsub(":|-","_",combp.male.annotated.top_10$region)
x = append(list("combp.male.annotated.top_10" = combp.male.annotated.top_10),tab)


writexl::write_xlsx(
  x = append(list("combp.male.annotated.top_10" = combp.male.annotated.top_10),tab),
  path = file.path("Michale/male_dist750/cnew.regions-p.bed_annotated_top_10.xlsx")
)


dim(combp.male)
dim(combp.male.annotated)
dim(combp.male.annotated.top_10)



combp.female.annotated <- combp.female %>% add_annotation
combp.female.annotated <- combp.female.annotated %>% mutate(
  n_probes = str_count(cpgs_in_region,"cg") 
)

writexl::write_xlsx(
  x = combp.female.annotated,
  path = file.path("Michale/female_750/cnew.regions-p.bed_annotated.xlsx")
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



# Add hyper hypo

library(readxl)
female_DMR_annotated <- read_excel("Michale/female_750/cnew.regions-p.bed_annotated.xlsx")

meta_df_AD_vs_CN.annotated.female <- readr::read_csv(
  file = paste0(
    dir.result.meta.analysis.ad,
    "FEMALE_meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg_annotated.csv"
  )
) %>% as.data.frame()

directions.female <- plyr::ldply(female_DMR_annotated$cpgs_in_region,.fun = function(x){
  
  
  estimate.sign.single <- meta_df_AD_vs_CN.annotated.female$direction.bacon[
    grep(gsub(",","|",x),meta_df_AD_vs_CN.annotated.female$cpg)
  ] %>% stringr::str_split(string = .,pattern = "") %>% unlist()
  percent.single.hypo <- sum(estimate.sign.single == "-")/length(estimate.sign.single)
  
  
  estimate.sign.meta <- ifelse(sign(meta_df_AD_vs_CN.annotated.female$estimate.bacon[
    grep(gsub(",","|",x),meta_df_AD_vs_CN.annotated.female$cpg)
  ]) == 1, "+","-")
  
  percent.meta.hypo <- sum(estimate.sign.meta == "-")/length(estimate.sign.meta)
  
  data.frame(
    "Direction_Meta_percent_estimate_hypo" = percent.meta.hypo,
    "Direction_Meta_estimate" =  paste(estimate.sign.meta,collapse = ""),
    "Direction_Single_cohort_percent_estimate_hypo" =  percent.single.hypo,
    "Direction_Single_cohort_estimate" =  paste(estimate.sign.single,collapse = "")
  )
})


female_DMR_annotated_with_direction <- cbind(female_DMR_annotated,directions.female) 
writexl::write_xlsx(female_DMR_annotated_with_direction,"Michale/female_750/cnew.regions-p.bed_annotated_with_directions.xlsx")
female_DMR_annotated_with_same_direction <- female_DMR_annotated_with_direction %>%
  filter(Direction_Meta_percent_estimate_hypo %in% c(0,1) & Direction_Single_cohort_percent_estimate_hypo %in% c(0,1))
writexl::write_xlsx(female_DMR_annotated_with_same_direction,"Michale/female_750/cnew.regions-p.bed_annotated_with_directions_final_filtered.xlsx")
female_DMR_annotated_with_same_direction <- female_DMR_annotated_with_same_direction %>%
  filter(z_p < 1E-5)
writexl::write_xlsx(female_DMR_annotated_with_same_direction,"Michale/female_750/cnew.regions-p.bed_annotated_with_directions_final_z_p_filtered.xlsx")


male_DMR_annotated <- read_excel("Michale/male_750/cnew.regions-p.bed_annotated.xlsx")

meta_df_AD_vs_CN.annotated.male <- readr::read_csv(
  file = paste0(
    dir.result.meta.analysis.ad,
    "MALE_meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg_annotated.csv"
  )
) %>% as.data.frame()



directions.male <- plyr::ldply(male_DMR_annotated$cpgs_in_region,.fun = function(x){
  
  
  estimate.sign.single <- meta_df_AD_vs_CN.annotated.male$direction.bacon[
    grep(gsub(",","|",x),meta_df_AD_vs_CN.annotated.male$cpg)
  ] %>% stringr::str_split(string = .,pattern = "") %>% unlist()
  percent.single.hypo <- sum(estimate.sign.single == "-")/length(estimate.sign.single)
  
  
  estimate.sign.meta <- ifelse(sign(meta_df_AD_vs_CN.annotated.male$estimate.bacon[
    grep(gsub(",","|",x),meta_df_AD_vs_CN.annotated.male$cpg)
  ]) == 1, "+","-")
  
  percent.meta.hypo <- sum(estimate.sign.meta == "-")/length(estimate.sign.meta)
  
  data.frame(
    "Direction_Meta_percent_estimate_hypo" = percent.meta.hypo,
    "Direction_Meta_estimate" =  paste(estimate.sign.meta,collapse = ""),
    "Direction_Single_cohort_percent_estimate_hypo" =  percent.single.hypo,
    "Direction_Single_cohort_estimate" =  paste(estimate.sign.single,collapse = "")
  )
})

male_DMR_annotated_with_direction <- cbind(male_DMR_annotated,directions.male)
writexl::write_xlsx(male_DMR_annotated_with_direction,"Michale/male_750/cnew.regions-p.bed_annotated_with_directions.xlsx")
male_DMR_annotated_with_same_direction <- male_DMR_annotated_with_direction %>%
  filter(Direction_Meta_percent_estimate_hypo %in% c(0,1) & Direction_Single_cohort_percent_estimate_hypo %in% c(0,1))
writexl::write_xlsx(male_DMR_annotated_with_same_direction,"Michale/male_750/cnew.regions-p.bed_annotated_with_directions_final_filtered.xlsx")
male_DMR_annotated_with_same_direction <- male_DMR_annotated_with_same_direction %>%
  filter(z_p < 1E-5)
writexl::write_xlsx(male_DMR_annotated_with_same_direction,"Michale/male_750/cnew.regions-p.bed_annotated_with_directions_final_z_p_filtered.xlsx")


writexl::write_xlsx(
  list(
    "Males DMRs" =  
      male_DMR_annotated_with_same_direction %>% dplyr::select(
        c("region","Direction_Meta_estimate","Direction_Single_cohort_estimate", "n_probes","z_p","z_sidak_p","GREAT_annotation","UCSC_RefGene_Name")
      ) %>% dplyr::arrange(z_sidak_p),
    "Females DMRs" =  female_DMR_annotated_with_same_direction %>% dplyr::select(
      c("region","Direction_Meta_estimate","Direction_Single_cohort_estimate","n_probes","z_p","z_sidak_p","GREAT_annotation","UCSC_RefGene_Name")
    ) %>% dplyr::arrange(z_sidak_p)
  ),path = "analysis_results/comb-p/table3_withSmokeProbes_z_p_filtered.xlsx"
)


# (3) make a venn diagram of overlapping filtered DMRs from males and females

male_DMR_annotated_with_same_direction <- readxl::read_xlsx(
  "~/TBL Dropbox/Wei Zhang//AD-meta-analysis-blood-samples-bySex/Michale/male_750/cnew.regions-p.bed_annotated_with_directions_final_z_p_filtered.xlsx"
)
female_DMR_annotated_with_same_direction <- readxl::read_xlsx(
  "~/TBL Dropbox/Wei Zhang//AD-meta-analysis-blood-samples-bySex/Michale/female_750/cnew.regions-p.bed_annotated_with_directions_final_z_p_filtered.xlsx"
)

gplots::venn(list("Comp-b Males" = male_DMR_annotated_with_same_direction$region,"Comp-b females" = female_DMR_annotated_with_same_direction$region))

findOverlaps(
  male_DMR_annotated_with_same_direction$region %>% MethReg::make_granges_from_names(),
  female_DMR_annotated_with_same_direction$region %>% MethReg::make_granges_from_names()
)

findOverlaps(
  male_DMR_annotated_with_same_direction$region %>% MethReg::make_granges_from_names(),
  female_DMR_annotated_with_same_direction$region %>% MethReg::make_granges_from_names()
)

# (4) venn diagram for overlapping DMRs with single CpGs in file
sig.cpgs.females <- readxl::read_xlsx("DRAFT-FIGURES-TABLES-REVISION_7-1-2022/Supp Table 2 sig_female_male_cpgs_AD_vs_CN_updated_final_LWedit - FINAL.xlsx",skip = 3,n_max = 24)
sig.cpgs.females <- sig.cpgs.females[-1,]
sig.cpgs.males <- readxl::read_xlsx("DRAFT-FIGURES-TABLES-REVISION_7-1-2022/Supp Table 2 sig_female_male_cpgs_AD_vs_CN_updated_final_LWedit - FINAL.xlsx", skip = 29)
colnames(sig.cpgs.males) <- colnames(sig.cpgs.females)
sig.cpgs.females <- sig.cpgs.females$cpg
sig.cpgs.males <- sig.cpgs.males$cpg


# (5) venn diagram for previous DMR results and current DMR results

previous_male_DMR <- readxl::read_xlsx("analysis_results/comb-p/table3.xlsx", sheet = 1)
previous_female_DMR <- readxl::read_xlsx("analysis_results/comb-p/table3.xlsx", sheet = 2)

#ggpubr::ggarrange(
  #plotlist = list(
   # ggvenn::ggvenn(
   #   list(
   #     "Comp-b analysis \n DMR sig. in males \n with same direction" = male_DMR_annotated_with_same_direction$region,
   #     "Comp-b analysis \n DMR sig. in females \n with same direction" = female_DMR_annotated_with_same_direction$region),
   #   fill_color = c("#0073C2FF","#CD534CFF"),
   #   set_name_size = 3
   # ),
    ggpubr::ggarrange(
      plotlist = list(
        
        ggvenn::ggvenn(
          list(
            "Males \n CpGs within Sig. DMR in\n comb-p analysis" = intersect(meta_df_AD_vs_CN.annotated.male$cpg,male_DMR_annotated_with_same_direction$cpgs_in_region %>% stringr::str_split(",") %>% unlist %>% unique),
            "Males \n Sig CpGs \nin single cpg analysis" = sig.cpgs.males
          ),
          fill_color = c("#EFC000FF","#868686FF"),
          set_name_size = 3
        ),
        
        ggvenn::ggvenn(
          list(
            "Females \nCpGs within Sig. DMR in\n comb-p analysis" = intersect(meta_df_AD_vs_CN.annotated.female$cpg,female_DMR_annotated_with_same_direction$cpgs_in_region %>% stringr::str_split(",") %>% unlist %>% unique),
            "Females \nSig CpGs \nin single cpg analysis" = sig.cpgs.females
          ),
          fill_color = c("#EFC000FF","#868686FF"),
          
          set_name_size = 3
        )
      ),ncol = 2, nrow = 1, labels = "AUTO"
    )
 # )
#)
    
ggsave(
  filename = "~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/analysis_results/others/Venn_and_Scatterplot_male_vs_female_results.pdf",
  width = 12,height = 8
)

ggpubr::ggarrange(
  plotlist = list(
    ggvenn::ggvenn(
      list(
        "Males \n CpGs within Sig. DMR in\n comb-p analysis \n with same direction" = male_DMR_annotated_with_same_direction$region,
        "Males \n CpGs within Sig. DMR in\n comb-p analysis \n previous results" = previous_male_DMR$region),
      fill_color = c("#0073C2FF","#CD534CFF"),
      set_name_size = 3
    ),

    ggvenn::ggvenn(
      list(
        "Females \n CpGs within Sig. DMR in\n comb-p analysis \n with same direction" = female_DMR_annotated_with_same_direction$region,
        "Females \n CpGs within Sig. DMR in\n comb-p analysis \n previous results" = previous_female_DMR$region
      ),
      fill_color = c("#EFC000FF","#868686FF"),
      set_name_size = 3
    )
  ), ncol = 2, nrow = 1
)


meta_df <- left_join(meta_df_AD_vs_CN.annotated.female, meta_df_AD_vs_CN.annotated.male, by = "cpg")

ggpubr::ggarrange(
  plotlist = list(
   
    ggvenn::ggvenn(
      list(
        "Male CpGs" = sig.cpgs.males,
        "Female CpGs" = sig.cpgs.females
      ),
      set_name_size = 4,
      set_name_color = c("#0073C2FF","#CD534CFF"),
      fill_color = c("#0073C2FF","#CD534CFF")
    ),
    
    ggvenn::ggvenn(
      list(
        "Male DMRs" = male_DMR_annotated_with_same_direction$region,
        "Female DMRs" = female_DMR_annotated_with_same_direction$region),
      set_name_size = 4,
      set_name_color = c("#0073C2FF","#CD534CFF"),
      fill_color = c("#0073C2FF","#CD534CFF")
    ), 

    ggpubr::ggscatter(
      meta_df, x = "estimate.bacon.y", y = "estimate.bacon.x",
      size = 0.5, add = "reg.line", add.params = list(color = "blue"),
      xlab = "estimate from Male analysis", 
      ylab = "estimate from Female analysis", 
      title = "Sig CpG Meta-analysis",
      cor.coef = T,
      cor.coeff.args = list(method = "spearman", aes(label = ..r.label..),label.x.npc = 0.9)
    ),
    
    ggplot(meta_df, aes(x = -log10(pVal.final.bacon.y), y = -log10(pVal.final.bacon.x))) + 
      geom_point(size = 0.5) + 
      geom_smooth(method='lm', color = "blue") + 
      ggpubr::stat_cor(method = "spearman", aes(label = ..r.label..),label.x.npc = 0.9) +
      theme_classic() + xlab("-log10(pValues) from Male analysis") + 
      ylab("-log10(pValues) from Female analysis") + ggtitle("Sig CpG Meta-analysis") +
      theme_pubr()
    
  ), ncol = 2, nrow = 2, labels = "AUTO"
)


ggsave(
  filename = "~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/analysis_results/others/Venn_and_Scatterplot_male_vs_female_results.pdf",
  width = 12,height = 8
)

meta_df <- meta_df %>% filter(
  pVal.final.bacon.x < 1E-05 | pVal.final.bacon.y < 1E-05
)

l = rbind(c(1,NA,2),
          c(NA,3,NA))
l = rbind(c(1,1), c(3,3))
layout.show(l)

select_grobs <- function(lay) {
  id <- unique(c(t(lay))) 
  id[!is.na(id)]
} 


  g1 = ggvenn::ggvenn(
    list(
      "Male CpGs" = sig.cpgs.males,
      "Female CpGs" = sig.cpgs.females
    ),
    set_name_size = 4,
    set_name_color = c("#0073C2FF","#CD534CFF"),
    fill_color = c("#0073C2FF","#CD534CFF")
  )
  
  g2 = ggvenn::ggvenn(
    list(
      "Male DMRs" = male_DMR_annotated_with_same_direction$region,
      "Female DMRs" = female_DMR_annotated_with_same_direction$region),
    set_name_size = 4,
    set_name_color = c("#0073C2FF","#CD534CFF"),
    fill_color = c("#0073C2FF","#CD534CFF")
  )


g3 = ggpubr::ggscatter(
  meta_df, x = "estimate.bacon.x", y = "estimate.bacon.y",
  size = 1,
  xlab = "estimate from Female analysis", 
  ylab = "estimate from Male analysis", 
  title = "Sig CpG (pValue < 1E-05) Meta-analysis",
  cor.coef = T,
  cor.coeff.args = list(method = "spearman", aes(label = ..r.label..),label.x.npc = 0.9)
) 

(g1 +   g2) / (plot_spacer() +  g3 + theme(plot.margin = unit(c(.3,0,2,0), "inch")) + plot_spacer() )


ggsave(
  filename = "~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/analysis_results/others/Venn_cpgs_male_vs_female_results_A.pdf",
  plot = g1, 
  width = 6,height = 4
)
ggsave(
  filename = "~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/analysis_results/others/Venn_dmr_male_vs_female_results_B.pdf",
  plot = g2, 
  width = 6,height = 4
)
ggsave(
  filename = "~/TBL Dropbox/Wei Zhang/AD-meta-analysis-blood-samples-bySex/analysis_results/others/scatter_cpgs_male_vs_female_results_C.pdf",
  plot = g3, 
  width = 6,height = 4
)
