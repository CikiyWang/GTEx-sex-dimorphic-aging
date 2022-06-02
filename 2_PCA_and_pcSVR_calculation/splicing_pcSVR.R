# calculate AS pcSVR

# load packages
library("data.table")
require("parallel")

# define functions
euclidean_dist <- function(x,y){
  z <- as.numeric(dist(rbind(x,y),method = "euclidean"))
  return(z)
}
pcSVR <- function(X,Y){
  # row = gene; col=sample
  nF <- ncol(X)
  nM <- ncol(Y)
  X_bar <- rowMeans(X)
  Y_bar <- rowMeans(Y)
  Signal <- euclidean_dist(X_bar,Y_bar)
  
  sum_square_Fs <- 0
  for (i in 1:nF) {
    sum_square_F <- euclidean_dist(X[,i],X_bar)^2
    sum_square_Fs <- sum_square_Fs + sum_square_F
  }
  noise_F <- sum_square_Fs/(nF-1)
  
  sum_square_Ms <- 0
  for (i in 1:nM) {
    sum_square_M <- euclidean_dist(Y[,i],Y_bar)^2
    sum_square_Ms <- sum_square_Ms + sum_square_M
  }
  noise_M <- sum_square_Ms/(nM-1)
  
  pcSVR <- Signal/sqrt((noise_F/nF)+(noise_M/nM))
  return(pcSVR)
}
mean_norm <- function(x){
  mean_est <- mean(x)
  sd <- mean(abs(x-mean_est))
  z <- (x-mean_est)/sd
  return(z)
}
rowmin <- function(x){
  return(as.numeric(apply(x, 1, min)))
}
perm_pval <- function(x){
  pval <- 1-(length(x[x>1])/length(x))
  return(pval)
}
fastread <- function(file,sep){
  df <- fread(file,sep = sep)
  df <- as.data.frame(df)
  rownames(df) <- df$V1
  df <- df[,-1]
  return(df)
}

# directory
PSI_dir <- "" 
pheno_dir <- ""
imp_PSI_dir <- "" 
info_dir <- ""
out_dir <- ""

setwd(out_dir)

age_cutoff_low <- 40
age_cutoff_high <- 60
n.cores <- 4

tissue_dir <- ""
print(tissue_dir)

impute_PSI_filter <- fastread(paste(imp_PSI_dir,"/GTEx_v8_",tissue_dir,"_ase_psi_filter_tpm_pc_imp.txt",sep = ""),sep = "\t")

# filter samples
if(file.info(paste(pheno_dir,"/",tissue_dir,"_pheno.txt",sep = ""))$size <= 2){next}
pheno <- read.table(paste(pheno_dir,"/",tissue_dir,"_pheno.txt",sep = ""),sep = "\t")
dim(pheno) 
rownames(pheno) <- pheno$V1
colnames_PSI <- paste(unlist(lapply(strsplit(colnames(impute_PSI_filter),"-",fixed = T),function(x){x[1]})),
                      unlist(lapply(strsplit(colnames(impute_PSI_filter),"-",fixed = T),function(x){x[2]})),sep = "-")
intersect_samples <- intersect(pheno$V1,colnames_PSI)
pheno <- pheno[intersect_samples,]
colnames(pheno) <- c("ID","SEX","AGE","RACE","BMI","DTHHRDY")
pheno$age_class <- factor(ifelse(pheno$AGE<=age_cutoff_low,"1:Young",
                                 ifelse(pheno$AGE>age_cutoff_high,"3:Old","2:Middle")))
impute_PSI_filter <- impute_PSI_filter[,which(colnames_PSI %in% intersect_samples)]

# filter sex-chromosomal genes
info <- read.delim(paste(info_dir,"/GTEx_v8_",tissue_dir,"_ase_info.txt",sep = ""),header = T)
impute_PSI_filter <- impute_PSI_filter[!(rownames(impute_PSI_filter) %in% info[info$chr %in% c("chrX","chrY"),]$events_id),]


PSI_pca <- prcomp(t(impute_PSI_filter),scale. = T)

# calculate PCs which can explain 80% variation
sds <- c()
for (i in 1:length(PSI_pca$sdev)) {
  sds <- c(sds,sum((PSI_pca$sdev/sum(PSI_pca$sdev))[1:i]))
}
sel_pc <- length(sds[sds<0.8])+1


PSI_PC1 <- as.matrix(PSI_pca$x[,1:sel_pc])
colnames(PSI_PC1) <- paste("PC",1:sel_pc,sep = "")
PSI_PC1 <- t(PSI_PC1)
colnames_PSI_PC1 <- paste(unlist(lapply(strsplit(colnames(PSI_PC1),"-",fixed = T),function(x){x[1]})),
                          unlist(lapply(strsplit(colnames(PSI_PC1),"-",fixed = T),function(x){x[2]})),sep = "-")

PC_M <- as.matrix(PSI_PC1[,colnames_PSI_PC1 %in% rownames(pheno[pheno$SEX==1,])])
PC_F <- as.matrix(PSI_PC1[,colnames_PSI_PC1 %in% rownames(pheno[pheno$SEX==2,])])
PC_Sex <- cbind(PC_M,PC_F)
dim(PC_M)
dim(PC_F)
PC_O <- as.matrix(PSI_PC1[,colnames_PSI_PC1 %in% rownames(pheno[pheno$AGE>age_cutoff_high,])])
PC_Y <- as.matrix(PSI_PC1[,colnames_PSI_PC1 %in% rownames(pheno[pheno$AGE<=age_cutoff_low,])])
PC_Age <- cbind(PC_O,PC_Y)
dim(PC_Y)
dim(PC_O)

n_permutation <- 10000

# filter tissues with small sample sizes
print(c(ncol(PC_F),ncol(PC_M),ncol(PC_Y),ncol(PC_O)))
if ((ncol(PC_F)<=2) | (ncol(PC_M)<=2) | (ncol(PC_Y)<=2) | (ncol(PC_O)<=2)){next;}

# calculate pcSVRs with group lables
pcSVRs_tissue <- mclapply(1:n_permutation,function(x){
  PC_F_sample <- PC_F[,sample(1:ncol(PC_F),0.5*min(ncol(PC_F),ncol(PC_M)))]
  PC_M_sample <- PC_M[,sample(1:ncol(PC_M),0.5*min(ncol(PC_F),ncol(PC_M)))]
  pcSVR_sex_tissue <- pcSVR(PC_F_sample,PC_M_sample)
  
  PC_Y_sample <- PC_Y[,sample(1:ncol(PC_Y),0.5*min(ncol(PC_Y),ncol(PC_O)))]
  PC_O_sample <- PC_O[,sample(1:ncol(PC_O),0.5*min(ncol(PC_Y),ncol(PC_O)))]
  pcSVR_age_tissue <-pcSVR(PC_Y_sample,PC_O_sample)
  
  return(list(pcSVR_sex_tissue,pcSVR_age_tissue))
},mc.cores = n.cores)

pcSVRs_sex_tissue <- unlist(lapply(pcSVRs_tissue, function(x) return(x[[1]])))
pcSVRs_age_tissue <- unlist(lapply(pcSVRs_tissue, function(x) return(x[[2]])))

pcSVRs_sex_tissue_mean <- mean(pcSVRs_sex_tissue)
pcSVRs_age_tissue_mean <- mean(pcSVRs_age_tissue)

# calculate pcSVRs without group lables
pcSVRs_tissue <- mclapply(1:n_permutation,function(x){
  PC_1_sample <- PC_Sex[,sample(1:ncol(PC_Sex),0.5*min(ncol(PC_F),ncol(PC_M)))]
  PC_2_sample <- PC_Sex[,sample(1:ncol(PC_Sex),0.5*min(ncol(PC_F),ncol(PC_M)))]
  pcSVR_1_tissue <- pcSVR(PC_1_sample,PC_2_sample)
  
  PC_1_sample <- PC_Age[,sample(1:ncol(PC_Age),0.5*min(ncol(PC_Y),ncol(PC_O)))]
  PC_2_sample <- PC_Age[,sample(1:ncol(PC_Age),0.5*min(ncol(PC_Y),ncol(PC_O)))]
  pcSVR_2_tissue <-pcSVR(PC_1_sample,PC_2_sample)
  
  return(list(pcSVR_1_tissue,pcSVR_2_tissue))
},mc.cores = n.cores)

pcSVRs_1_tissue <- unlist(lapply(pcSVRs_tissue, function(x) return(x[[1]])))
pcSVRs_2_tissue <- unlist(lapply(pcSVRs_tissue, function(x) return(x[[2]])))

pcSVR_sex_perm_pval <- length(pcSVRs_1_tissue[pcSVRs_1_tissue>pcSVRs_sex_tissue_mean])/10000
pcSVR_age_perm_pval <- length(pcSVRs_2_tissue[pcSVRs_2_tissue>pcSVRs_age_tissue_mean])/10000

pcSVR_sex_pval <- perm_pval(pcSVRs_sex_tissue)
pcSVR_age_pval <- perm_pval(pcSVR_age_perm_pval)

print(c(pcSVRs_sex_tissue_mean,pcSVRs_age_tissue_mean))
print(c(pcSVR_sex_perm_pval,pcSVR_age_perm_pval))
