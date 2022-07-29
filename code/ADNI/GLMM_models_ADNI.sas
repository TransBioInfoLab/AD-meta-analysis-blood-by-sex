options mprint; 

%let dir = C:\DATA_for_SAS_models; 

%let res = C:\RESULTS; 

%macro p; 

%do i = 1 %to 100; 

PROC IMPORT OUT=ad_cn DATAFILE= "&dir\ADNI_1000_cpgs_with_smok_probes&i..csv" REPLACE; GETNAMES=yes; RUN;

data ad_cn; set ad_cn; beta_pct = beta * 100; run;

proc sort data = ad_cn; by cpg rid PlateNumber; run;

ods exclude all; ods noresults;

proc glimmix data = ad_cn; 
by cpg; 
	class RID PlateNumber; 

	model DIAGNOSIS (descending)  = beta_pct age_at_visit  PlateNumber B NK CD4T Mono Granulocytes / dist = binary s; 

	random int /subject = RID; 

	*ods output Tests3 = res_beta; 

	ods output ParameterEstimates = res_beta; 
run; 
ods exclude none; ods results;

data res_beta; set res_beta; if Effect = "beta_pct"; format Probt e10.; drop PlateNumber; run;

proc export data=res_beta dbms=csv outfile="&res\females_glmm_using_beta_&i..csv" replace; run;

%end; 

%mend; 

%p; 
