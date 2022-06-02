# Differential GE analysis

# load packages
library("SmartSVA")
library("edgeR")
library("dplyr")
require("parallel")
library("purrr")
library("data.table")

# define function
diff_exp_glm <- function(d,annot){
  
  colnames_d <- paste(unlist(lapply(strsplit(colnames(d),"-",fixed = T),function(x){x[1]})),
                      unlist(lapply(strsplit(colnames(d),"-",fixed = T),function(x){x[2]})),sep = "-")
  sample <- pheno[pheno$ID %in% colnames_d,]
  sample <- sample[intersect(colnames_d,sample$ID),c("SEX","age_class")]
  sample$SEX <- factor(sample$SEX)
  d <- d[,colnames_d %in% intersect(colnames_d,rownames(sample))]
  colnames_d <- paste(unlist(lapply(strsplit(colnames(d),"-",fixed = T),function(x){x[1]})),
                      unlist(lapply(strsplit(colnames(d),"-",fixed = T),function(x){x[2]})),sep = "-")
  sample <- sample[colnames_d,]
  
  d_norm　<- as.data.frame(t(apply(d, 1,mean_norm)))
  
  cat("\nStart calcuate Surrogate Variables ...\n")
  cat("======================================\n")
  
  Y.r <- t(resid(lm(t(d_norm) ~ SEX+age_class+SEX*age_class, data=sample)))
  n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
  mod <- model.matrix(~SEX+age_class+SEX*age_class, sample)
  sv.obj <- sva(as.matrix(d_norm), mod, mod0=NULL,method="irw",n.sv=n.sv)
  mod1 <- model.matrix(~SEX+age_class+SEX*age_class+sv.obj$sv,sample)
  
  cat("\n\nStart linear regression in 20 threads ...\n")
  cat("======================================\n")
  
  sv_lm_result <- mclapply(as.list(as.data.frame(t(d))),function(x){
    TPM <- data.frame(TPM=mean_norm(as.numeric(x)),SEX=factor(sample[,1]),age_class=sample[,2])
    df_for_lm <- cbind(TPM,sv.obj$sv)
    fit <- glm(TPM~.+SEX*age_class,data = df_for_lm)
    coef <- as.data.frame(coef(summary(fit)))
    coef_Estimate <- coef$Estimate
    coef_pval <- coef$`Pr(>|t|)`
    
    TPMs <- as.numeric(x)
    fc_TPM_sex <- mean(TPMs[which(df_for_lm$SEX==2)])/mean(TPMs[which(df_for_lm$SEX==1)])
    fc_TPM_age_middle <- mean(TPMs[which(df_for_lm$age_class=="2:Middle")])/mean(TPMs[which(df_for_lm$age_class=="1:Young")])
    fc_TPM_age_old <- mean(TPMs[which(df_for_lm$age_class=="3:Old")])/mean(TPMs[which(df_for_lm$age_class=="1:Young")])
    
    return(list(coef_Estimate,c(fc_TPM_sex,fc_TPM_age_old),coef_pval))
  },mc.cores = n.cores)
  
  coef_Estimates <- t(as.data.frame(lapply(sv_lm_result, function(x) return(x[[1]]))))
  foldchanges <- t(as.data.frame(lapply(sv_lm_result, function(x) return(x[[2]]))))
  coef_pvals <- t(as.data.frame(lapply(sv_lm_result, function(x) return(x[[3]]))))

  colnames(coef_Estimates) <- paste("coef",colnames(mod1),sep="_")
  colnames(foldchanges) <- c("fc_TPM_sex","fc_TPM_age_old")
  colnames(coef_pvals) <- paste("pval",colnames(mod1),sep="_")
  rownames(coef_Estimates) <- rownames(d)
  rownames(foldchanges) <- rownames(d)
  rownames(coef_pvals) <- rownames(d)
  
  coef_Estimates <- coef_Estimates[,!(colnames(coef_Estimates) %in% c("coef_(Intercept)",paste("coef_sv.obj$sv",1:n.sv,sep = "")))]
  coef_pvals <- coef_pvals[,!(colnames(coef_pvals) %in% c("pval_(Intercept)",paste("pval_sv.obj$sv",1:n.sv,sep = "")))]
  
  diff_list <- as.data.frame(cbind(coef_Estimates,coef_pvals,foldchanges))
  diff_list$symbol <- annot[rownames(d),]$gene_name
  return(diff_list)
}
mean_norm <- function(x){
  mean_est <- mean(x)
  sd <- mean(abs(x-mean_est))
  z <- (x-mean_est)/sd
  return(z)
}
fastread <- function(file,sep){
  df <- fread(file,sep = sep)
  df <- as.data.frame(df)
  rownames(df) <- df$V1
  df <- df[,-1]
  return(df)
}

# directory
pheno_dir <- ""
TPM_dir <- "" 
TPM_info_dir <- "" 
out_dir <- ""
setwd(out_dir)

age_cutoff_low <- 40
age_cutoff_high <- 60
n.cores <- 4

tissue_dir <- ""
print(tissue_dir)

# TPM matrix
TPM <- fastread(paste(TPM_dir,"/GTEx_v8_",tissue_dir,"_exp_tpm.txt",sep = ""),sep = "\t")
dim(TPM)

info <- read.table(paste(TPM_info_dir,"/GTEx_v8_",tissue_dir,"_exp_tpm_info.txt",sep = ""),sep = "\t",header = T)
info_unique <- info %>% unique()
rownames(info_unique) <- info_unique$gene_id

# filter exp
TPM_count <- TPM[which(rowMeans(TPM)>1),]
dim(TPM_count)

pheno <- read.table(paste(pheno_dir,"/",tissue_dir,"_pheno.txt",sep = ""),sep = "\t")
dim(pheno)
colnames_TPM_count <- paste(unlist(lapply(strsplit(colnames(TPM_count),"-",fixed = T),function(x){x[1]})),
                            unlist(lapply(strsplit(colnames(TPM_count),"-",fixed = T),function(x){x[2]})),sep = "-")
pheno <- pheno[pheno$V1 %in% colnames_TPM_count,]
rownames(pheno) <- pheno$V1
colnames(pheno) <- c("ID","SEX","AGE","RACE","BMI","DTHHRDY")
pheno$age_class <- factor(ifelse(pheno$AGE<=age_cutoff_low,"1:Young",
                                 ifelse(pheno$AGE>age_cutoff_high,"3:Old","2:Middle")))

group_nrow <- pheno %>% 
  group_by(SEX,age_class) %>%
  do(data.frame(nrow=nrow(.)))

group_nrow_sex<- pheno %>% 
  group_by(SEX) %>%
  do(data.frame(nrow=nrow(.)))

group_nrow_age<- pheno %>% 
  group_by(age_class) %>%
  do(data.frame(nrow=nrow(.)))

if(length(group_nrow_sex$SEX)==1){next}
if(length(group_nrow_age$age_class)==1){next}
if(min(group_nrow_sex$nrow)<3){next}
if(min(group_nrow_age$nrow)<3){next}
if(ncol(TPM_count)<20){next}
if(length(levels(factor(ifelse(pheno$AGE<=age_cutoff_low,"1:Young",ifelse(pheno$AGE>age_cutoff_high,"3:Old","2:Middle")))))==1){next}

anno_tissue <- info_unique[rownames(TPM_count),]

diff_list <- diff_exp_glm(TPM_count,anno_tissue)

write.table(diff_list,paste(out_dir,"/GTEx_v8_",tissue_dir,"_diff_exp_inter_MN_padj.txt",sep=""),sep="\t",col.names=T,quote=F)

