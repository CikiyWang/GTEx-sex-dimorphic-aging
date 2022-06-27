######################################
# sBASEs for AD prediction
######################################

library("sda")
library("tm")
library("mixOmics")
library("dplyr")
library("mlbench")
library("caret")
library("caTools")

# annotation
pc_ENSG_symbol <- read.table("gencode.v36.protein_coding_gene_ENSG_symbol.bed",sep = "\t")

diff_splicing_dir <- ""
diff_sex_splicing_dir <- ""
info_dir <- ""

tissue_dir <- "Brain_DECISION"

# input diff results
diff_results <- read.table(paste(diff_splicing_dir,"/GTEx_v8_",tissue_dir,"_diff_splice_inter_padj.txt",sep=""),sep = "\t",row.names = 1,header = T)
interactions <- diff_results[(diff_results$pval_SEX2.age_class3.Old < 0.05),]
diff_sex_results <- read.table(paste(diff_sex_splicing_dir,"/GTEx_v8_",tissue_dir,"_diff_splice_inter_padj_sex_stratified.txt",sep=""),sep = "\t",header = T)
diff_splicing_female <- diff_sex_results[(diff_sex_results$pval_age_class3.Old_female<0.05) & (abs(diff_sex_results$deltapsi_PSI_age_old_female)>0.05),]
diff_splicing_male <- diff_sex_results[(diff_sex_results$pval_age_class3.Old_male<0.05) & (abs(diff_sex_results$deltapsi_PSI_age_old_male)>0.05),]
female_sp <- diff_sex_results[setdiff(rownames(diff_splicing_female),rownames(diff_splicing_male)),]
male_sp <- diff_sex_results[setdiff(rownames(diff_splicing_male),rownames(diff_splicing_female)),]
sb_female_sp <- diff_results[intersect(rownames(interactions),rownames(female_sp)),]
sb_male_sp <- diff_results[intersect(rownames(interactions),rownames(male_sp)),]
diff_AD_Ctrl <- read.table("GSE174367_FC_ase_psi_filter_tpm_pc_imp_diff_splicing_sex_stratified.txt",
                           sep = "\t",row.names = 1,header = T)
splicing_merge <- merge(diff_sex_results,diff_AD_Ctrl,by="row.names",suffixes = c("aging","AD"))
rownames(splicing_merge) <- splicing_merge$Row.names
splicing_merge <- splicing_merge[,-1]

# read matrix, pheno, and split to different sexes
PSI_dir <- ""
pheno_dir <- ""
GTEx_PSI <- read.table(paste(PSI_dir,"/GTEx_v8_",tissue_dir,"_ase_psi_filter_tpm_pc_imp.txt",sep = ""),sep = "\t",header = T,row.names = 1)
colnames(GTEx_PSI) <- gsub(".","-",colnames(GTEx_PSI),fixed = T)
dim(GTEx_PSI)
pheno <- read.table(paste(pheno_dir,"/",tissue_dir,"_pheno.txt",sep = ""),sep = "\t")
dim(pheno)
pheno <- pheno[pheno$V3>60 | pheno$V3<=40,]
rownames(pheno) <- pheno$V1
colnames_GTEx_PSI <- paste(unlist(lapply(strsplit(colnames(GTEx_PSI),"-",fixed = T),function(x){x[1]})),
                           unlist(lapply(strsplit(colnames(GTEx_PSI),"-",fixed = T),function(x){x[2]})),sep = "-")
intersect_samples <- intersect(pheno$V1,colnames_GTEx_PSI)
pheno <- pheno[intersect_samples,]
colnames(pheno) <- c("ID","SEX","AGE","RACE","BMI","DTHHRDY")
pheno$age_class <- factor(ifelse(pheno$AGE<=40,"1:Young",
                                 ifelse(pheno$AGE>60,"3:Old","2:Middle")))
GTEx_PSI <- GTEx_PSI[,which(colnames_GTEx_PSI %in% intersect_samples)]
colnames_GTEx_PSI <- paste(unlist(lapply(strsplit(colnames(GTEx_PSI),"-",fixed = T),function(x){x[1]})),
                           unlist(lapply(strsplit(colnames(GTEx_PSI),"-",fixed = T),function(x){x[2]})),sep = "-")
GTEx_PSI_M <- GTEx_PSI[,which(colnames_GTEx_PSI %in% rownames(pheno[pheno$SEX==1,]))]
GTEx_PSI_F <- GTEx_PSI[,which(colnames_GTEx_PSI %in% rownames(pheno[pheno$SEX==2,]))]
pheno_M <- pheno[paste(unlist(lapply(strsplit(colnames(GTEx_PSI_M),"-",fixed = T),function(x){x[1]})),
                       unlist(lapply(strsplit(colnames(GTEx_PSI_M),"-",fixed = T),function(x){x[2]})),sep = "-"),]
pheno_F <- pheno[paste(unlist(lapply(strsplit(colnames(GTEx_PSI_F),"-",fixed = T),function(x){x[1]})),
                       unlist(lapply(strsplit(colnames(GTEx_PSI_F),"-",fixed = T),function(x){x[2]})),sep = "-"),]

# read matrix in AD dataset 
AD_PSI <- read.table("GSE174367_FC_ase_psi_filter_tpm_pc_imp.txt",sep = "\t",header = T,row.names = 1)
meta <- read.csv("GSE174367_all_metadata.txt",sep = ",",header = T,row.names = 1)
colnames(meta)
sample <- meta[colnames(AD_PSI),c("AGE","diagnosis","sex")]
AD_PSI_F <- AD_PSI[,rownames(sample[sample$sex=="female",])]
AD_PSI_M <- AD_PSI[,rownames(sample[sample$sex=="male",])]
sample_F <- sample[colnames(AD_PSI_F),] 
sample_M <- sample[colnames(AD_PSI_M),]

GTEx_PSI_F_sel <- GTEx_PSI_F[rownames(splicing_merge),]
GTEx_PSI_M_sel <- GTEx_PSI_M[rownames(splicing_merge),]
AD_PSI_F_sel <- AD_PSI_F[rownames(splicing_merge),]
AD_PSI_M_sel <- AD_PSI_M[rownames(splicing_merge),]

GTEx_PSI_F_sBASEs <- GTEx_PSI_F[rownames(GTEx_PSI_F) %in% rownames(female_sp),]
GTEx_PSI_M_sBASEs <- GTEx_PSI_M[rownames(GTEx_PSI_M) %in% rownames(male_sp),]
AD_PSI_F_sBASEs <- AD_PSI_F[rownames(AD_PSI_F) %in% rownames(female_sp),]
AD_PSI_M_sBASEs <- AD_PSI_M[rownames(AD_PSI_M) %in% rownames(male_sp),]

AD_PSI_F_exsBASEs <- AD_PSI_F[!(rownames(AD_PSI_F) %in% rownames(female_sp)),]
AD_PSI_M_exsBASEs <- AD_PSI_M[!(rownames(AD_PSI_M) %in% rownames(male_sp)),]

# randomly selected AS events in sex-stratified model
AD_PSI_F_random <- AD_PSI_F_exsBASEs[sample(1:nrow(AD_PSI_F_exsBASEs),nrow(AD_PSI_F_sBASEs)),]
AD_PSI_M_random  <- AD_PSI_M_exsBASEs[sample(1:nrow(AD_PSI_M_exsBASEs),nrow(AD_PSI_M_sBASEs)),]

# sBASEs in merge-sexes model
AD_PSI_MBASEs <- AD_PSI[rownames(AD_PSI) %in% rownames(male_sp),]
AD_PSI_FBASEs <- AD_PSI[rownames(AD_PSI) %in% rownames(female_sp),]

###########################################################################################################################
############################################# sex-stratified model ########################################################
###########################################################################################################################

######################################
# start training in females
######################################
setwd("./female")
dataset = as.data.frame(t(AD_PSI_F_sBASEs[,]))
dataset$Y=factor(sample_F$diagnosis)

#set.seed(123)

df_split <- data.frame()
for (i in 1:100) {
  split = sample.split(dataset$Y, SplitRatio = 0.9)
  df_split <- rbind(df_split,split)
}
colnames(df_split) <- 1:ncol(df_split)
write.table(df_split,paste("Female_rf_cv_all_foldrep.txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = T)

df_performance <- data.frame()
selVars_lens <- c()
library(pROC)
for (i in 1:100) { 
  # outer loop
  print(i)
  training_set = dataset[df_split[i,] == TRUE,]
  test_set = dataset[df_split[i,] == FALSE,]
  
  # feature selection
  rfFuncs$summary <- twoClassSummary
  control <- rfeControl(functions=rfFuncs, method="cv", number=5)
  # run the RFE algorithm
  results <- rfe(training_set[,-ncol(training_set)], training_set$Y, sizes=seq(10,50,by=5), rfeControl=control,metric = "ROC")
  
  rf.pred.training=predict(results,test_set,type="prob")
  confussion.training<-confusionMatrix(factor(rf.pred.training$pred), test_set$Y)
  
  varimp_data <- data.frame(feature = row.names(varImp(results)),
                            importance = varImp(results,scale = FALSE)[,1])
  write.table(varimp_data,paste("Female_rf_cv_innerloop_varImp_foldrep_",i,".txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = T)

  selVars <- predictors(results)
  print(length(selVars))
  write.table(selVars,paste("Female_rf_cv_innerloop_selVars_foldrep_",i,".txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = T)
  write.table(results$results,paste("Female_rf_cv_innerloop_results_foldrep_",i,".txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = T)
  
  pROC_obj <- roc(test_set$Y,rf.pred.training$AD,
                  smoothed = F,direction=">",
                  ci=F, ci.alpha=0.5, stratified=FALSE,
                  plot=F, auc.polygon=TRUE,max.auc.polygon=T,grid=TRUE,auc.polygon.col="#E6BAA9",
                  print.auc=TRUE, show.thres=TRUE)
  
  auc <- pROC_obj$auc
  
  df_res <- as.data.frame(t(as.matrix(c(confussion.training$overall,confussion.training$byClass))))
  df_res$auc <- pROC_obj$auc

  df_performance <- rbind(df_performance,df_res)
}
write.table(df_performance,"Female_rf_ncv_performance.txt",sep = "\t",quote = F,row.names = T,col.names = T)


######################################
# start training in males
######################################
## male
setwd("./male")
dataset = as.data.frame(t(AD_PSI_M_sBASEs[,]))
dataset$Y=factor(sample_M$diagnosis)

df_split <- data.frame()
for (i in 1:100) {
  split = sample.split(dataset$Y, SplitRatio = 0.9)
  df_split <- rbind(df_split,split)
}
colnames(df_split) <- 1:ncol(df_split)
write.table(df_split,paste("Male_rf_cv_all_foldrep.txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = T)

df_performance <- data.frame()
selVars_lens <- c()
library(pROC)
for (i in 1:100) { 
  # outer loop
  print(i)
  training_set = dataset[df_split[i,] == TRUE,]
  test_set = dataset[df_split[i,] == FALSE,]
  
  # feature selection
  rfFuncs$summary <- twoClassSummary
  control <- rfeControl(functions=rfFuncs, method="cv",number = 5)
  trainctrl <- trainControl(classProbs= TRUE,
                            summaryFunction = twoClassSummary)
  # run the RFE algorithm
  results <- rfe(training_set[,-ncol(training_set)], training_set$Y, sizes=seq(10,50,by=5), rfeControl=control)
  
  rf.pred.training=predict(results,test_set,type="prob")
  confussion.training<-confusionMatrix(factor(rf.pred.training$pred), test_set$Y)
  
  varimp_data <- data.frame(feature = row.names(varImp(results,scale = FALSE)),
                            importance = varImp(results,scale = FALSE)[,1])
  write.table(varimp_data,paste("Male_rf_cv_innerloop_varImp_foldrep_",i,".txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = T)

  selVars <- predictors(results)
  print(length(selVars))
  write.table(selVars,paste("Male_rf_cv_innerloop_selVars_foldrep_",i,".txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = T)
  write.table(results$results,paste("Male_rf_cv_innerloop_results_foldrep_",i,".txt",sep = ""),sep = "\t",quote = F,row.names = T,col.names = T)
  
  pROC_obj <- roc(test_set$Y,rf.pred.training$AD,
                  smoothed = F,direction=">",
                  ci=F, ci.alpha=0.5, stratified=FALSE,
                  plot=F, auc.polygon=TRUE,max.auc.polygon=T,grid=TRUE,auc.polygon.col="#E6BAA9",
                  print.auc=TRUE, show.thres=TRUE)
  
  auc <- pROC_obj$auc
  
  df_res <- as.data.frame(t(as.matrix(c(confussion.training$overall,confussion.training$byClass))))
  df_res$auc <- pROC_obj$auc
  
  df_performance <- rbind(df_performance,df_res)
}
write.table(df_performance,"Male_rf_ncv_performance.txt",sep = "\t",quote = F,row.names = T,col.names = T)
