# filter psi values 
# 1. filter samples with > 50% NA events
# 2. remove na < 5%
# 3. JC > 10
# 4. constant 
# 5. max-min > 0.05
# 6. sd > 0.01
# 7. range in 0.05 ~ 0.95
# 8. filter tpm>1 and protein coding genes
# 9 imputation by KNN

library("data.table")
require("parallel")
library("impute")

rowmin <- function(x){
  return(as.numeric(apply(x, 1, min)))
}

fastread <- function(file,sep){
  df <- fread(file,sep = sep)
  df <- as.data.frame(df)
  rownames(df) <- df$V1
  df <- df[,-1]
  return(df)
}

PSI_dir <- "" 
IJC_dir <- "" 
SJC_dir <- "" 
CI_dir <- "" 
pheno_dir <- ""
imp_PSI_dir <- "" 
info_dir <- ""
TPM_dir <- "" 
Preprocessing_dir <- "" 

age_cutoff_low <- 40
age_cutoff_high <- 60
n.cores <- 10

tissue_dir <- ""
print(tissue_dir)

## read tables
# PSI matrix
PSI <- fastread(paste(PSI_dir,"/GTEx_v8_",tissue_dir,"_ase_psi.txt",sep = ""),
                sep = "\t")
dim(PSI)
# ijc matrix
IJC <- fastread(paste(IJC_dir,"/GTEx_v8_",tissue_dir,"_ase_ijc.txt",sep = ""),
                sep = "\t")
dim(IJC)
# sjc matrix
SJC <- fastread(paste(SJC_dir,"/GTEx_v8_",tissue_dir,"_ase_sjc.txt",sep = ""),
                sep = "\t")
dim(SJC)
JC <- IJC+SJC
# sjc matrix
CI <- fastread(paste(CI_dir,"/GTEx_v8_",tissue_dir,"_ase_ci.txt",sep = ""),
               sep = "\t")
dim(SJC)

# filter samples with phenotypes
if(file.info(paste(pheno_dir,"/",tissue_dir,"_pheno.txt",sep = ""))$size <= 2){next}
pheno <- read.table(paste(pheno_dir,"/",tissue_dir,"_pheno.txt",sep = ""),sep = "\t")
dim(pheno) 
rownames(pheno) <- pheno$V1
colnames_PSI <- paste(unlist(lapply(strsplit(colnames(PSI),"-",fixed = T),function(x){x[1]})),
                      unlist(lapply(strsplit(colnames(PSI),"-",fixed = T),function(x){x[2]})),sep = "-")
intersect_samples <- intersect(pheno$V1,colnames_PSI)
pheno <- pheno[intersect_samples,]
colnames(pheno) <- c("ID","SEX","AGE","RACE","BMI","DTHHRDY")
pheno$age_class <- factor(ifelse(pheno$AGE<=age_cutoff_low,"1:Young",
                                 ifelse(pheno$AGE>age_cutoff_high,"3:Old","2:Middle")))
PSI <- PSI[,which(colnames_PSI %in% intersect_samples)]
JC <- JC[,which(colnames_PSI %in% intersect_samples)]
print(dim(PSI))
PSI[PSI == -1] <- NA

## Start filtering
# filter samples with missing values
na_count_sample <- mclapply(as.data.frame(PSI),function(x){return(length(x[is.na(x)]))},mc.cores = n.cores)
data_na_sample_0.5 <- PSI[,as.numeric(which(unlist(na_count_sample)<(0.5*nrow(PSI))))]
dim(data_na_sample_0.5)
df_step1 <- as.data.frame(summary(factor(unlist(lapply(strsplit(rownames(data_na_sample_0.5),"@"),function(x){x[1]})))))

# filter events with missing values
na_count <- mclapply(as.data.frame(t(data_na_sample_0.5)),function(x){return(length(x[is.na(x)]))},mc.cores = n.cores)
data_na_0.05 <- data_na_sample_0.5[which(na_count<(0.05*ncol(data_na_sample_0.5))),]
df_step2 <- as.data.frame(summary(factor(unlist(lapply(strsplit(rownames(data_na_0.05),"@"),function(x){x[1]})))))

# filter events with JC > 20
new_JC <- JC[rownames(data_na_0.05),]
data_na_0.05 <- data_na_0.05[which(rowMeans(new_JC)>10),]
df_step3 <- as.data.frame(summary(factor(unlist(lapply(strsplit(rownames(data_na_0.05),"@"),function(x){x[1]})))))

if (nrow(data_na_0.05)==0){next;}
if (ncol(data_na_0.05)==0){next;}

# filter events constantly
levels <- mclapply(as.data.frame(t(data_na_0.05)),function(x){return(length(levels(factor(x))))},mc.cores = n.cores)
data_na_0.05_levels <- data_na_0.05[which(unlist(levels)>2),]
df_step4 <- as.data.frame(summary(factor(unlist(lapply(strsplit(rownames(data_na_0.05_levels),"@"),function(x){x[1]})))))

# filter max-min > 0.05
max_min <- apply(data_na_0.05_levels, 1, function(x){return(max(x[is.na(x)==F])-min(x[is.na(x)==F]))})
data_max_min_0.01 <- data_na_0.05_levels[which(max_min>0.05),]
df_step5 <- as.data.frame(summary(factor(unlist(lapply(strsplit(rownames(data_max_min_0.01),"@"),function(x){x[1]})))))

# filter sd > 0.01
sds <- apply(data_max_min_0.01, 1,function(x){return(sd(x[is.na(x)==F]))})
data_max_min_0.01_var <- data_max_min_0.01[rownames(data_max_min_0.01) %in% names(which(sds>0.01)),]
df_step6 <- as.data.frame(summary(factor(unlist(lapply(strsplit(rownames(data_max_min_0.01_var),"@"),function(x){x[1]})))))

# filter events in 0.05~0.95
data_max_min_0.01_var_range <- data_max_min_0.01_var[which(rowMeans(data_max_min_0.01_var,na.rm = T) < 0.95 & rowMeans(data_max_min_0.01_var,na.rm = T) > 0.05),]
df_step7 <- as.data.frame(summary(factor(unlist(lapply(strsplit(rownames(data_max_min_0.01_var_range),"@"),function(x){x[1]})))))
filter_psi <- data_max_min_0.01_var_range

# filter TPM and pc
TPM <- fastread(paste(TPM_dir,"/GTEx_v8_",tissue_dir,"_exp_tpm.txt",sep = ""),sep = "\t")
dim(TPM)
TPM_count <- TPM[which(rowMeans(TPM)>1),]
dim(TPM_count)

# filter protein coding genes
pc_ENSG_symbol <- read.table("/picb/rnasys2/wangsiqi/Annotation/human/GRCh38/gencode.v36.protein_coding_gene_ENSG_symbol.bed",sep = "\t")
tpm_pc_symbol <- pc_ENSG_symbol[pc_ENSG_symbol$V5 %in% rownames(TPM_count),]$V6
info <- read.delim(paste(info_dir,"/GTEx_v8_",tissue_dir,"_ase_info.txt",sep = ""),header = T)
anno_tpm <- info[which(info$symbol %in% tpm_pc_symbol),]
dim(info)
dim(anno_tpm)
anno_tissue <- anno_tpm[anno_tpm$full_id %in% rownames(filter_psi),]
filter_psi_tpm_pc <- filter_psi[anno_tissue$full_id,]
dim(anno_tissue)
dim(filter_psi_tpm_pc)
print(dim(filter_psi_tpm_pc))
print(summary(factor(unlist(lapply(strsplit(rownames(filter_psi_tpm_pc),"@"),function(x){x[1]})))))
df_step8 <- as.data.frame(summary(factor(unlist(lapply(strsplit(rownames(filter_psi_tpm_pc),"@"),function(x){x[1]})))))

# inputation by KNN
imputation=impute.knn(as.matrix(filter_psi_tpm_pc), k=10, maxp=dim(filter_psi_tpm_pc)[1])
impute_PSI_filter <- imputation$data
impute_PSI_filter[which(impute_PSI_filter>1)]=1
impute_PSI_filter[which(impute_PSI_filter<0)]=0
print(dim(impute_PSI_filter))

# output PSI matrix
write.table(impute_PSI_filter,paste(imp_PSI_dir,"/GTEx_v8_",tissue_dir,"_ase_psi_filter_tpm_pc_imp.txt",sep = ""),
            sep = "\t",row.names = T,col.names = T,quote = F)

# stat numbers of AS events in each step
df_ase_type_count <- cbind(df_step1,df_step2,df_step3,df_step4,
                           df_step5,df_step6,df_step7,df_step8)
colnames(df_ase_type_count) <- c("Raw","NA_5%",
                                 "Mean_JC_10","Constant_Value",
                                 "Max_min_0.05","sd_0.01",
                                 "Mean_0.05_0.95","Pc_TPM")
write.table(df_ase_type_count,paste(Preprocessing_dir,"/GTEx_v8_",tissue_dir,"_ase_type_stat.txt",sep = ""),
            sep = "\t",row.names = T,col.names = T,quote = F)
