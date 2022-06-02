# Differential AS analysis

# load packages
library("SmartSVA")
library("edgeR")
library("purrr")
library("data.table")
library("dplyr")
require("parallel")

# define function
diff_ase_glm <- function(d,annot){
  
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
    PSI <- data.frame(PSI=mean_norm(as.numeric(x)),SEX=factor(sample[,1]),age_class=factor(sample[,2]))
    df_for_lm <- cbind(PSI,sv.obj$sv)
    fit <- glm(PSI~.+SEX*age_class,data = df_for_lm,na.action = na.omit)
    coef <- as.data.frame(coef(summary(fit)))

    coef_Estimate <- coef$Estimate
    coef_pval <- coef$`Pr(>|t|)`
    
    PSIs <- as.numeric(x)
    deltapsi_PSI_sex <- mean(PSIs[which(df_for_lm$SEX==2)], na.rm = TRUE)-mean(PSIs[which(df_for_lm$SEX==1)], na.rm = TRUE)
    deltapsi_PSI_age_old <- mean(PSIs[which(df_for_lm$age_class=="3:Old")], na.rm = TRUE)-mean(PSIs[which(df_for_lm$age_class=="1:Young")], na.rm = TRUE)
    
    return(list(coef_Estimate,c(deltapsi_PSI_sex,deltapsi_PSI_age_old),coef_pval))
  },mc.cores = n.cores)
  
  coef_Estimates <- t(as.data.frame(lapply(sv_lm_result, function(x) return(x[[1]]))))
  deltapsi <- t(as.data.frame(lapply(sv_lm_result, function(x) return(x[[2]]))))
  coef_pvals <- t(as.data.frame(lapply(sv_lm_result, function(x) return(x[[3]]))))
  
  colnames(coef_Estimates) <- paste("coef",colnames(mod1),sep="_")
  colnames(deltapsi) <- c("deltapsi_PSI_sex","deltapsi_PSI_age_old")
  colnames(coef_pvals) <- paste("pval",colnames(mod1),sep="_")
  rownames(coef_Estimates) <- rownames(d)
  rownames(deltapsi) <- rownames(d)
  rownames(coef_pvals) <- rownames(d)
  
  coef_Estimates <- coef_Estimates[,!(colnames(coef_Estimates) %in% c("coef_(Intercept)",paste("coef_sv.obj$sv",1:n.sv,sep = "")))]
  coef_pvals <- coef_pvals[,!(colnames(coef_pvals) %in% c("pval_(Intercept)",paste("pval_sv.obj$sv",1:n.sv,sep = "")))]
  
  diff_list <- as.data.frame(cbind(coef_Estimates,coef_pvals,deltapsi))
  diff_list$symbol <- annot$geneSymbol
  diff_list$full_id <- annot$full_id
 
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
rowmin <- function(x){
  return(as.numeric(apply(x, 1, min)))
}

pheno_dir <- ""
imp_PSI_dir <- "" 
info_dir <- ""
out_dir <- ""

age_cutoff_low <- 40
age_cutoff_high <- 60
n.cores <- 2

setwd(out_dir)

tissue_dir <- ""
print(tissue_dir)

# PSI matrix
imp_psi <- fastread(paste(imp_PSI_dir,"/GTEx_v8_",tissue_dir,"_ase_psi_filter_tpm_pc_imp.txt",sep = ""),sep="\t")

if(file.info(paste(pheno_dir,"/",tissue_dir,"_pheno.txt",sep = ""))$size <= 2){next}
pheno <- read.table(paste(pheno_dir,"/",tissue_dir,"_pheno.txt",sep = ""),sep = "\t")
dim(pheno) 
rownames(pheno) <- pheno$V1
colnames_PSI <- paste(unlist(lapply(strsplit(colnames(imp_psi),"-",fixed = T),function(x){x[1]})),
                      unlist(lapply(strsplit(colnames(imp_psi),"-",fixed = T),function(x){x[2]})),sep = "-")
intersect_samples <- intersect(pheno$V1,colnames_PSI)
pheno <- pheno[intersect_samples,]
colnames(pheno) <- c("ID","SEX","AGE","RACE","BMI","DTHHRDY")
pheno$age_class <- factor(ifelse(pheno$AGE<=age_cutoff_low,"1:Young",
                                 ifelse(pheno$AGE>age_cutoff_high,"3:Old","2:Middle")))
imp_psi <- imp_psi[,which(colnames_PSI %in% intersect_samples)]
print(dim(imp_psi))

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
if(ncol(imp_psi)<20){next}
if(length(levels(factor(ifelse(pheno$AGE<=age_cutoff_low,"1:Young",ifelse(pheno$AGE>age_cutoff_high,"3:Old","2:Middle")))))==1){next}

anno <- read.delim(paste(info_dir,"/GTEx_v8_",tissue_dir,"_ase_info.txt",sep = ""),header = T)
anno$geneSymbol <- anno$symbol
anno$event_class <- anno$events_type
anno$chr <- anno$chr
rownames(anno) <- anno$full_id
anno_tissue <- anno[rownames(imp_psi),]
dim(imp_psi)
dim(anno)
dim(anno_tissue)

diff_list <- diff_ase_glm(imp_psi,anno_tissue)

write.table(diff_list,paste(out_dir,"/GTEx_v8_",tissue_dir,"_diff_splice_inter_padj.txt",sep=""),sep="\t",col.names=T,quote=F)



