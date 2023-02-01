# GTEx-sex-dimorphic-aging

the source code of GTEx sex dimorphic aging project

R version 4.0.4 (2021-02-15)

### Required R packages:
#### R packages:
	'impute=1.64.0';
	'parallel=4.0.4'
	'dplyr=1.2.2'
	'limma=3.46.0'
	'edgeR=3.32.1'
	'sva=3.38.0'
	'caret=6.0-93'
	'caTools=1.18.2'
	'zoo=1.8-10'
	'forecast=8.16'
#### Functional enrichment analysis was performed by R packages:
	'msigdbr=7.5.1'
	'clusterProfiler=3.18.1'
	'DOSE=3.16.0'.

### Usage: 
	Rscript *.R or in RStudio


## Data preprocessing:
**Processing the quantification of gene expression (GE) and alternative splicing (AS) as judged by Paean**

>**Gene expression:**  
>>TPM > 1 and protein coding genes  

>**Alternative splicing:**  
>>**(1)** averaged total reads counts in spliced junctions (i.e., spliced in counts + spliced out counts) > 10  
  **(2)** non-constant PSI values across samples  
  **(3)** max(PSI) – min(PSI) > 0.05  
  **(4)** standard deviation > 0.01  
  **(5)** averaged PSI value in the range from 0.05 to 0.95  
  **(6)** TPM of the spliced genes >1  

>**input:** data matrix with genes/AS events in rows and sample name in columns  
**output:** processed data matrix  


## PCA and pcSVR calculation:
**Calculation of PCA and PCA-based Signal-to-Variation Ratio (pcSVR) to evaluate the contributions to global transcriptomic variations**

>**input:** GE/AS matrix and phenotype file (sex, age, etc.)  
**output:** print pcSVR value and corresponding empirical p-value  
 

## Differential analysis:
**sex-/age-differential analysis and sex-stratified differential analysis**

>**input:** GE/AS matrix and phenotype file (sex, age, etc.)  
**output:** results of differential analysis in .txt file 

>**columns:** 
>>coef_SEX2, coef_age_class2:Middle, coef_age_class3:Old, coef_SEX2:age_class2:Middle, coef_SEX2:age_class3:Old
pval_age_class2:Middle, pval_age_class3:Old, pval_SEX2:age_class2:Middle, pval_SEX2:age_class3:Old
fc_TPM_sex, fc_TPM_age_old, symbol


## AD association analysis:
**RandomForest classifier for predicting Alzheimer's Disease (AD) with feature selections procedure (RFE)**

>**Sex-stratified model:**  
  >>sBASEs in females predict female AD; sBASEs in males predicte male AD;  

>**Control model:**  
  >>Sex-stratified model by randomly selected AS events  
  Merge-sexes model: sBASEs in females/males predict AD  

>**input:** AS matrix with sBASEs in rows and phenotype file (e.g., AD, Control)  
**output:** Performance of AD prediction in 100 iterations  


## SFs-AS events Network:
**Construction the regulatory networks between splicing factors (SFs) and AS events**

>**Step1:** Spearman's correlation between SFs expression (TPM) and AS events (PSI) < 0.05  
**Step2:** Significant regulation by SF Knockdown (shRNA-seq from ENCODE)  
**Step3:** Existence of the SFs motif signals in the adjacent region (+/-300nt) of splice sites  

>**Step1:** Step1_correlation_SFs_and_AS_evnets.R: 
>>input GE and AS data matrix; 
>>output spearman's rho and p-value;  

>**Step2:** Step2_diff_shRBP_AS_events_ENCODE.R: 
>>input PSI matrix of shRNA-seq; 
>>output: differential AS events between control and RBP knockdown;  

>**Step3:** Step3_motif_identification_make_seqs_for_deepbind.R: 
>>input list of AS events and genome annotation; 
>>output the subsequences (sliding window=40nt) in the range from upstream to downstream 300nt of each AS event to run deepbind;  

>**Step4:** integration of the results from the 3 steps above:  
>>integration_networks_of_3_steps.R  


## Time series and Breakpoint analysis:
**Identification of chronological genes and AS events during aging**  
**Calculation breakpoints during aging**  
**Identification of Aging-determiant Genes (ADGs) which can significantly alter the aging breakpoints**  

>**Time_series_analysis_for_chronological.R:** 
>>input GE/AS matrix; output chronological GE/AS during aging.  

>**breakpoint_analysis_subsampled_genes_for_ADGs.R:** 
>>input hronological GE/AS during aging; output breakpoint and fitted TPM/PSI values.  
