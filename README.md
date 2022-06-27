# GTEx-sex-dimorphic-aging

the source code of GTEx sex dimorphic aging project

**1_Data_preprocessing:**
Processing the quantification of gene expression (GE) and alternative splicing (AS) as judged by Paean
GE: 
  TPM > 1 and protein coding genes
AS: 
  (1) averaged total reads counts in spliced junctions (i.e., spliced in counts + spliced out counts) > 10
  (2) non-constant PSI values across samples
  (3) max(PSI) – min(PSI) > 0.05
  (4) standard deviation > 0.01
  (5) averaged PSI value in the range from 0.05 to 0.95
  (6) TPM of the spliced genes >1

**2_PCA_and_pcSVR_calculation**
Calculation of PCA and PCA-based Signal-to-Variation Ratio (pcSVR) to evaluate the contributions to global transcriptomic variations

**3_Differential_analysis**
sex-/age-differential analysis and sex-stratified differential analysis

**4_AD_association**
RandomForest classifier for predicting Alzheimer's Disease (AD) with feature selections procedure (RFE)
Sex-stratified model:
  sBASEs in females predict female AD; sBASEs in males predicte male AD;
Control model:
  Sex-stratified model by randomly selected AS events
  Merge-sexes model: sBASEs in females/males predict AD

**5_SFs_RNA_Network_Construction**
Construction the regulatory networks between splicing factors (SFs) and AS events
Step1: Spearman's correlation between SFs expression (TPM) and AS events (PSI) < 0.05
Step2: Significant regulation by SF Knockdown (shRNA-seq from ENCODE)
Step3: Existence of the SFs motif signals in the adjacent region (+/-300nt) of splice sites
Integration code

**6_Time_series_and_Breakpoint_analysis**
Identification of chronological genes and AS events during aging
Calculation breakpoints during aging
Identification of Aging-determiant Genes (ADGs) which can significantly alter the aging breakpoints

