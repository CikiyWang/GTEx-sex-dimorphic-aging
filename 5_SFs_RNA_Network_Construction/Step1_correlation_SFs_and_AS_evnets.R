# Splicing factors regulation 
# step1
# spearman correlation between differential SF expression and AS events

######################################
# loading packages: 
######################################
library("ggsci")
library("ggpubr")
library("pheatmap")
library("dplyr")
require("parallel")

fastread <- function(file,sep){
  df <- fread(file,sep = sep)
  df <- as.data.frame(df)
  rownames(df) <- df$V1
  df <- df[,-1]
  return(df)
}

######################################
# directory: 
######################################
pheno_dir <- ""
imp_PSI_dir <- "" 
info_dir <- ""
TPM_dir <- "" 
diff_splice_dir <- ""
diff_exp_dir <- ""
out_cor_dir <- ""

out_dir <- ""
setwd(out_dir)

tissue_dir <- ""
print(tissue_dir)

######################################
# read data: 
######################################
SF <- read.table("GO_RNA_SPLICING_SFs.txt",sep = "\t",header = F)$V1
meta <- read.delim("experiment_report_full_metadata.tsv",sep = "\t",header = T)
shRBPs <- unique(meta$Target.of.assay)
focus_SFs <- intersect(SF,shRBPs)
pc_ENSG_symbol <- read.table("gencode.v36.protein_coding_gene_ENSG_symbol.bed",sep = "\t")
rownames(pc_ENSG_symbol) <- pc_ENSG_symbol$V5

diff_splicing_results <- read.table(paste(diff_splice_dir,"/GTEx_v8_",tissue_dir,"_diff_splice_inter_padj_sex_stratified.txt",sep=""),sep = "\t",row.names = 1,header = T)
diff_splicing_age_female <- diff_splicing_results[(diff_splicing_results$pval_age_class3.Old_female < 0.05) & (abs(diff_splicing_results$deltapsi_PSI_age_old_female) > 0.05),]
diff_splicing_age_male <- diff_splicing_results[(diff_splicing_results$pval_age_class3.Old_male < 0.05) & (abs(diff_splicing_results$deltapsi_PSI_age_old_male) > 0.05),]

# PSI
PSI <- fastread(paste(imp_PSI_dir,"/GTEx_v8_",tissue_dir,"_ase_psi_filter_tpm_pc_imp.txt",sep = ""),"\t")
dim(PSI)
# TPM
TPM <- fastread(paste(TPM_dir,"/GTEx_v8_",tissue_dir,"_exp_tpm.txt",sep = ""),"\t")
TPM <- TPM[rowMeans(TPM)>1,]
dim(TPM)
# phenotype
pheno <- fastread(paste(pheno_dir,"/",tissue_dir,"_pheno.txt",sep = ""),"\t")
dim(pheno)

samples <- intersect(colnames(PSI),colnames(TPM))
samples_id <- paste(unlist(lapply(strsplit(samples,"-",fixed = T),function(x){x[1]})),
                    unlist(lapply(strsplit(samples,"-",fixed = T),function(x){x[2]})),sep = "-")
male_samples <- samples[samples_id %in% rownames(pheno[pheno$V2==1,])]
female_samples <- samples[samples_id %in% rownames(pheno[pheno$V2==2,])]

# feamle 
PSI_diff_female <- PSI[diff_splicing_age_female$full_id_female,female_samples]
TPM_diff_female <- TPM[rownames(TPM) %in% pc_ENSG_symbol[pc_ENSG_symbol$V6 %in% focus_SFs,]$V5,female_samples]
symbols_female <- pc_ENSG_symbol[rownames(TPM_diff_female),]$V6
dim(PSI_diff_female)
dim(TPM_diff_female)

# male
PSI_diff_male <- PSI[diff_splicing_age_male$full_id_male,male_samples]
TPM_diff_male <- TPM[rownames(TPM) %in% pc_ENSG_symbol[pc_ENSG_symbol$V6 %in% focus_SFs,]$V5,male_samples]
symbols_male <- pc_ENSG_symbol[rownames(TPM_diff_male),]$V6
dim(PSI_diff_male)
dim(TPM_diff_male)


for (choose_sex in c("female","male")) {
  print(choose_sex)
  if(choose_sex == "female"){
    PSI_diff <- PSI_diff_female
    TPM_diff <- TPM_diff_female
    symbols <- symbols_female
  }
  if(choose_sex == "male"){
    PSI_diff <- PSI_diff_male
    TPM_diff <- TPM_diff_male
    symbols <- symbols_male
  }
  PSI_diff_var <- apply(PSI_diff, 1, var)
  TPM_diff_var <- apply(TPM_diff, 1, var)
  PSI_diff <- PSI_diff[PSI_diff_var>0,]
  TPM_diff <- TPM_diff[TPM_diff_var>0,]
  
  if(nrow(PSI_diff)==0 || nrow(TPM_diff)==0){next;}
  
  ######################################################################
  # step1 : filter correlation
  ######################################################################
  # correlation plot
  cor_res <- mclapply(as.data.frame(t(TPM_diff)),function(x){
    cors_r <- c()
    cors_pval <- c()
    for (n in 1:nrow(PSI_diff)) {
      res <- cor.test(as.numeric(x),as.numeric(PSI_diff[n,]),method = "spearman")
      cors_r <- c(cors_r,res$estimate)
      cors_pval <- c(cors_pval,res$p.value)
      #print(r)
    }
    return(list(as.numeric(cors_r),as.numeric(cors_pval)))
  },mc.cores = 4)
  
  df_cors_r <- t(as.data.frame(lapply(cor_res, function(x) return(x[[1]]))))
  df_cors_pval <- t(as.data.frame(lapply(cor_res, function(x) return(x[[2]]))))
  colnames(df_cors_r) <- rownames(PSI_diff)
  colnames(df_cors_pval) <- rownames(PSI_diff)
  rownames(df_cors_r) <- symbols
  rownames(df_cors_pval) <- symbols
  
  write.table(df_cors_r,paste(out_cor_dir,"/GTEx_v8_",tissue_dir,"_cor_rho_",choose_sex,".txt",sep = ""),sep = "\t",quote = F,col.names = T,row.names = T)
  write.table(df_cors_pval,paste(out_cor_dir,"/GTEx_v8_",tissue_dir,"_cor_pval_",choose_sex,".txt",sep = ""),sep = "\t",quote = F,col.names = T,row.names = T)
  
}


