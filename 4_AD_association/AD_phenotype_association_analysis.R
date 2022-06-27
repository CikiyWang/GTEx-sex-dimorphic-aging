# AD phenotype association analysis: tangles and plaque
require("parallel")
setwd("./")

metadata <- read.table("GSE174367_all_metadata.txt",sep = ",",header = T)
sample_id <- read.table("GSM_id_Total_bulk_RNAseq.tsv",sep = "\t",header = T)
sample_clinical <- read.delim("Bulk_RNA-seq_metadata.tsv",sep = "\t",header = T)

TPM <- read.table("GSE174367_FC_exp_tpm.txt",sep = "\t",row.names = 1,header = T)
TPM <- TPM[rowMeans(TPM)>1,]
imp_PSI <- read.table("GSE174367_FC_ase_psi_filter_tpm_pc_imp.txt",
                      sep = "\t",header = T)
meta_FC <- metadata[metadata$Run %in% colnames(imp_PSI),]
sample_id_FC <- sample_id[sample_id$id %in% meta_FC$Sample.Name,]
sample_clinical$SampleID <- paste("Sample",sample_clinical$Sample.ID,sep = "-")
sample_clinical_FC <- sample_clinical[sample_clinical$SampleID %in% unique(sample_id_FC$sample),]
rownames(sample_clinical_FC) <- sample_clinical_FC$SampleID
rownames(meta_FC) <- meta_FC$Sample.Name
rownames(sample_id_FC) <- sample_id_FC$id
sample_info <- merge(meta_FC,sample_id_FC,by="row.names")
rownames(sample_info) <- sample_info$sample
sample_info <- sample_info[,-1]
pheno <- merge(sample_info,sample_clinical_FC,by="row.names")[,c("Run","id","sample","Age","Sex","Tangle.Stage",
                                                                 "Plaque.Stage","Diagnosis","RIN","Sequencing.Group")]
rownames(pheno) <- pheno$Run
pheno$Tangle.Stage[is.na(pheno$Tangle.Stage)] <- NA
pheno$Tangle.Stage <- as.numeric(gsub("Stage ","",pheno$Tangle.Stage))

pheno$Plaque.Stage[pheno$Plaque.Stage=="None"] <- 0
pheno$Plaque.Stage[pheno$Plaque.Stage=="Stage A"] <- 1
pheno$Plaque.Stage[pheno$Plaque.Stage=="Stage B"] <- 2
pheno$Plaque.Stage[pheno$Plaque.Stage=="Stage C"] <- 3
pheno$Plaque.Stage <- as.numeric(pheno$Plaque.Stage)

pheno <- pheno[colnames(imp_PSI),]
pheno_M <- pheno[pheno$Sex == "M",]
pheno_F <- pheno[pheno$Sex == "F",]

imp_PSI_F <- imp_PSI[,rownames(pheno_F)]
imp_PSI_M <- imp_PSI[,rownames(pheno_M)]
dim(imp_PSI_F)
dim(imp_PSI_M)


mean_norm <- function(x){
  mean_est <- mean(x)
  sd <- mean(abs(x-mean_est))
  z <- (x-mean_est)/sd
  return(z)
}

tangles_association <- function(d){
  d_res <- mclapply(as.data.frame(t(d)), function(x){
    sample <- pheno[colnames(d),]
    df_for_lm <- sample[,c("Age","Tangle.Stage")]
    df_for_lm$PSI <- mean_norm(as.numeric(x))
    
    fit <- lm(PSI~.,data = df_for_lm)
    coef <- as.data.frame(coef(summary(fit)))
    
    coef_Estimate <- coef$Estimate
    coef_pval <- coef$`Pr(>|t|)`
    return(list(coef_Estimate,coef_pval))
  },mc.cores = 1)
  coef_Estimates <- t(as.data.frame(lapply(d_res, function(x) return(x[[1]]))))
  coef_pvals <- t(as.data.frame(lapply(d_res, function(x) return(x[[2]]))))
  tangles_pvals <- as.numeric(coef_pvals[,3])
  tangles_Estimates <- as.numeric(coef_Estimates[,3])
  df <- data.frame(events=rownames(d),pvals=tangles_pvals,Estimates=tangles_Estimates)
  return(df)
}
imp_PSI_F_tangles <- tangles_association(imp_PSI_F)
imp_PSI_M_tangles <- tangles_association(imp_PSI_M)

plaque_association <- function(d){
  d_res <- mclapply(as.data.frame(t(d)), function(x){
    sample <- pheno[colnames(d),]
    df_for_lm <- sample[,c("Age","Plaque.Stage")]
    df_for_lm$PSI <- mean_norm(as.numeric(x))
    
    fit <- lm(PSI~.,data = df_for_lm)
    coef <- as.data.frame(coef(summary(fit)))
    
    coef_Estimate <- coef$Estimate
    coef_pval <- coef$`Pr(>|t|)`
    return(list(coef_Estimate,coef_pval))
  },mc.cores = 1)
  coef_Estimates <- t(as.data.frame(lapply(d_res, function(x) return(x[[1]]))))
  coef_pvals <- t(as.data.frame(lapply(d_res, function(x) return(x[[2]]))))
  plaque_pvals <- as.numeric(coef_pvals[,3])
  plaque_Estimates <- as.numeric(coef_Estimates[,3])
  df <- data.frame(events=rownames(d),pvals=plaque_pvals,Estimates=plaque_Estimates)
  return(df)
}
imp_PSI_F_plaque <- plaque_association(imp_PSI_F)
imp_PSI_M_plaque <- plaque_association(imp_PSI_M)

write.table(imp_PSI_F_tangles,"tangles_F_AS_events_pvals.txt",quote = F,sep = "\t",row.names = F,col.names = T)
write.table(imp_PSI_M_tangles,"tangles_M_AS_events_pvals.txt",quote = F,sep = "\t",row.names = F,col.names = T)
write.table(imp_PSI_F_plaque,"plaque_F_AS_events_pvals.txt",quote = F,sep = "\t",row.names = F,col.names = T)
write.table(imp_PSI_M_plaque,"plaque_M_AS_events_pvals.txt",quote = F,sep = "\t",row.names = F,col.names = T)

