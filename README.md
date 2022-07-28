# Integrative analyses revealed sex-specific DNA methylation differences that underlie Alzheimer’s disease

Tiago C. Silva, Wei Zhang, Juan I. Young,  Lissette Gomez, Michael A. Schmidt,  Xi Chen, Eden R. Martin, Lily Wang

# Citing this repository
[![DOI](https://zenodo.org/badge/457491267.svg)](https://zenodo.org/badge/latestdoi/457491267)

### Description

This github repository includes scripts used for the analyses in the above manuscript. 

In this work, we performed a sex-specific meta-analysis of two large independent blood-based epigenome-wide association studies, the ADNI and AIBL studies (see references below), with a total of 1284 whole blood samples (633 female samples and 651 male samples). Furthermore, to identify blood-based DNA methylation markers that also change with underlying neuropathology in the brain, we next performed a sex-specific cross-tissue meta-analysis by combining these blood DNA methylation datasets with four additional DNA methylation datasets, which included a total of  1,030 prefrontal cortex brain samples (642 female samples and 388 male samples). Moreover, we used two complementary analytical strategies within each dataset, a sex-stratified analysis that examined methylation to AD associations in male and female samples separately, and a methylation-by-sex interaction analysis that compared the magnitude of these associations between different sexes. Our findings highlighted distinct sex-specific epigenetic architectures underlie AD and provide a useful resource for future biomarker studies in AD.  

### 1. Study cohorts, Preprocessing of DNA methylation data, Single Cohort analysis


| File                 | Dataset | Link |
|----------------------|-------------|-------------|
| ADNI/ADNI_SAS_bySex.Rmd  |   ADNI  | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood-by-sex/blob/main/code/ADNI/ADNI_SAS_bySex.Rmd) |
| AIBL/AIBL_bySex.Rmd           |   AIBL    | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood-by-sex/blob/main/code/AIBL/AIBL_bySex.rmd) |
| Matched_data_ADNI/matched_RNA_DNAm_data_and_residuals_bySex.R          |   ADNI    | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood-by-sex/blob/main/code/Matched_data_ADNI/matched_RNA_DNAm_data_and_residuals_bySex.R) |
| Clinical/clinical_info.Rmd          |   ADNI & AIBL & addneuromed   | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood-by-sex/blob/main/code/Clinical/clinical_info.Rmd) |
| Clinical/clinical_info_brain.Rmd          |  GASPARONI & LONDON & MtSinai & ROSMAP  | [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood-by-sex/blob/main/code/Clinical/clinical_info_brain.Rmd) |


### 2. Blood samples meta-analysis

| File                 | Link |
|----------------------|-------------|
| meta-analysis/meta-analysis-two-cohorts_glm_by_sex.Rmd      |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood-by-sex/blob/main/code/meta-analysis/meta-analysis-two-cohorts_glm_by_sex.Rmd) |
| meta-analysis-two-cohorts-interaction-glm.Rmd      |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood-by-sex/blob/main/code/meta-analysis/meta-analysis-two-cohorts-interaction-glm.Rmd) |

#### 2.1 Blood samples meta-analysis results 

| Result | File               | Link |
|---------|--------------------|-------------|
|Female| MALE_meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg.csv     |  [Link](https://github.com/TransBioInfoLab/AD-meta-analysis-blood-by-sex/blob/main/results/MALE_meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg.csv) |
|Male| FEMALE_meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg.csv     |  [Link](https://github.com/TransBioInfoLab/AD-meta-analysis-blood-by-sex/blob/main/results/FEMALE_meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_single_cpg.csv) |
|Interaction| meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_interaction_single_cpg.csv     |  [Link](https://github.com/TransBioInfoLab/AD-meta-analysis-blood-by-sex/blob/main/results/meta_analysis_glm_fixed_effect_ADNI_and_AIBL_AD_vs_CN_interaction_single_cpg.csv) |


### 3. Cross-tissue meta-analysis	

| File                 | Link |
|----------------------|-------------|
| cross_tissue_meta_analysis/cross_tissue_meta_analysis_male.Rmd        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood-by-sex/blob/main/code/cross_tissue_meta_analysis/cross_tissue_meta_analysis_male.Rmd) |
| cross_tissue_meta_analysis/cross_tissue_meta_analysis_female.Rmd        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood-by-sex/blob/main/code/cross_tissue_meta_analysis/cross_tissue_meta_analysis_female.Rmd) |

#### 3.1 Cross-tissue meta-analysis results

| Result | File               | Link |
|---------|--------------------|-------------|
|Male| MALE_cross_tissue_meta_analysis_glm_using_AD_vs_CN_single_cpg.csv   |  [Link](https://github.com/TransBioInfoLab/AD-meta-analysis-blood-by-sex/blob/main/results/MALE_cross_tissue_meta_analysis_glm_using_AD_vs_CN_single_cpg.csv) |

### 4. Correlations between methylation levels of significant CpGs and DMRs in AD with expressions of nearby genes

| File                 | Link |
|----------------------|-------------|
| DNAm_vs_RNA/Blood_ADNI_RNA_vs_cpg_bySex.R      |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood-by-sex/blob/main/code/DNAm_vs_RNA/Blood_ADNI_RNA_vs_cpg_bySex.R) |
| DNAm_vs_RNA/Brain_ROSMAP_RNA_vs_cpg_bySex.R      |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood-by-sex/blob/main/code/DNAm_vs_RNA/Brain_ROSMAP_RNA_vs_cpg_bySex.R) |

### 5. Out-of-sample validations of AD-associated DNAm differences in an external cohort - Methylation_Risk_scores

| File                 | Link |
|----------------------|-------------|
| Methylation_Risk_scores/Methylation_risk_scores_both.Rmd        |  [Link to the script](https://github.com/TransBioInfoLab/AD-meta-analysis-blood-by-sex/blob/main/code/Methylation_Risk_scores/Methylation_risk_scores_both.Rmd ) |


# For reproducible research

The following R packages are required: 

```r
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install(version = "3.14",ask = FALSE) # Install last version of Bioconductor

list.of.packages <- c(
  "bacon",
  "EpiSmokEr",
  "DMRcate",                                      
  "doParallel",                                   
  "dplyr",                                        
  "DT",                                           
  "EpiDISH",                                      
  "ExperimentHub",                                
  "fgsea",                                        
  "GenomicRanges",                                
  "GEOquery",                                     
  "ggpubr",                                       
  "ggrepel",                                      
  "gridExtra",                                    
  "gt",                                           
  "GWASTools",                                    
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "lubridate",                                    
  "lumi",                                         
  "meta",                                         
  "metap",                                        
  "MethReg",                                      
  "minfi",                                        
  "missMethyl",                                   
  "mygene",                                       
  "plyr",                                         
  "readr",                                        
  "readxl",                                       
  "ReMapEnrich",                                  
  "RPMM",                                         
  "RVenn",                                        
  "sm",                                           
  "stats",                                        
  "SummarizedExperiment",                         
  "tidyr",                                        
  "wateRmelon",                                   
  "writexl" 
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

devtools::install_github("igordot/msigdbr")
```

For ADNIMERGE, download it from https://ida.loni.usc.edu/: Merged ADNI 1/GO/2 Packages for R

```r
install.packages("/path/to/ADNIMERGE_0.0.1.tar.gz", repos = NULL, type = "source")
```
The platform information are: 

```r
version  R version 4.2.0 (2022-04-22)
 os       macOS Big Sur 11.4          
 system   x86_64, darwin17.0          
 ui       RStudio                     
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
 tz       America/New_York            
 date     2021-07-12      
```

# Acknowledgement
Data used in preparation of this article were obtained from the Alzheimer’s Disease Neuroimaging Initiative
(ADNI) database (adni.loni.usc.edu). As such, the investigators within the ADNI contributed to the design
and implementation of ADNI and/or provided data but did not participate in analysis or writing of this report. A complete listing of ADNI investigators can be found at:
http://adni.loni.usc.edu/wp-content/uploads/how_to_apply/ADNI_Acknowledgement_List.pdf

# References

1. Vasanthakumar, A. et al. Harnessing peripheral DNA methylation differences in the Alzheimer's Disease Neuroimaging Initiative (ADNI) to reveal novel biomarkers of disease. Clin Epigenetics 12, 84 (2020).

2. Ellis, K.A. et al. Enabling a multidisciplinary approach to the study of ageing and Alzheimer's disease: an update from the Australian Imaging Biomarkers and Lifestyle (AIBL) study. Int Rev Psychiatry 25, 699-710 (2013).
