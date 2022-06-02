######################################################################
# step2 : filter diff AS events with RBP binding signals by deepbind
######################################################################

##################################################
# load packages:
##################################################
library("SmartSVA")
library("edgeR")
library("purrr")
library("data.table")
library("dplyr")
require("parallel")

##################################################
# define functions:
##################################################
fastread <- function(file,sep){
  df <- fread(file,sep = sep)
  df <- as.data.frame(df)
  rownames(df) <- df$V1
  df <- df[,-1]
  return(df)
}
diff_ase_glm <- function(d,annot){
  
  # sva
  sample <- data.frame(group = factor(c("G1","G1","G2","G2")))
  rownames(sample) <- colnames(d)
  
  cat("\nStart calcuate Surrogate Variables ...\n")
  cat("======================================\n")
  d <- na.omit(d)
  num <- apply(d, 1, function(x){return(length(levels(factor(x))))})
  d <- d[num>1,]
  rownames(annot) <- anno$event_class
  annot <- annot[rownames(d),]
  
  mod <- model.matrix(~., sample)
  
  cat("\n\nStart linear regression in 20 threads ...\n")
  cat("======================================\n")
  sv_lm_result <- mclapply(as.list(as.data.frame(t(d))),function(x){
    PSI <- data.frame(PSI=as.numeric(x),group = factor(c("G1","G1","G2","G2")))
    df_for_lm <- PSI
    fit <- glm(PSI~.,data = df_for_lm)
    coef <- as.data.frame(coef(summary(fit)))
    coef$fdr <- p.adjust(coef$`Pr(>|t|)`,method = "fdr")
    
    coef_Estimate <- coef$Estimate
    coef_fdr <- coef$fdr
    
    PSIs <- as.numeric(x)
    deltapsi_PSI <- mean(PSIs[which(df_for_lm$group=="G2")])-mean(PSIs[which(df_for_lm$group=="G1")])
    
    return(list(coef_Estimate,coef_fdr,c(deltapsi_PSI)))
  },mc.cores = n.cores)
  
  coef_Estimates <- t(as.data.frame(lapply(sv_lm_result, function(x) return(x[[1]]))))
  deltapsi <- t(as.data.frame(lapply(sv_lm_result, function(x) return(x[[2]]))))
  
  colnames(coef_Estimates) <- paste("coef",colnames(mod),sep="_")
  colnames(deltapsi) <- c("deltapsi_PSI")
  rownames(coef_Estimates) <- rownames(d)
  rownames(deltapsi) <- rownames(d)
  
  coef_Estimates <- coef_Estimates[,!(colnames(coef_Estimates) %in% c("coef_(Intercept)"))]
  coef_Estimates <- as.data.frame(coef_Estimates)
  colnames(coef_Estimates) <- "coef_groupG2"

  diff_list <- as.data.frame(cbind(coef_Estimates,coef_fdrs,deltapsi))

  diff_list$symbol <- annot$geneSymbol
  diff_list$events <- annot$event_class
  
  return(diff_list)
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

##################################################
# parameters:
##################################################
PSI_dir <- ""
IJC_dir <- ""
SJC_dir <- ""
info_dir <- ""
out_dir <- ""

setwd(out_dir)
RBP_ids <- gsub("_ase_psi.txt","",list.files(PSI_dir,pattern = "*_ase_psi.txt"))
n.cores <- 10

for (RBP_id in RBP_ids) {
  print(RBP_id)

  imp_psi <- fastread(paste(PSI_dir,"/",RBP_id,"_ase_psi.txt",sep = ""),"\t")
  anno <- read.delim(paste(info_dir,"/",RBP_id,"_ase_info.txt",sep = ""),header = T)

  anno$geneSymbol <- unlist(lapply(strsplit(anno$ID,":"),function(x){x[length(x)-1]}))
  anno$event_class <- unlist(lapply(strsplit(anno$ID,":"),function(x){x[length(x)]}))
  anno$chr <- unlist(lapply(strsplit(anno$ID,":"),function(x){x[1]}))
  
  diff_list <- diff_ase_glm(imp_psi,anno)
  
  write.table(diff_list,paste(out_dir,"/splicing/",RBP_id,"_diff_splice_inter_padj.txt",sep=""),sep="\t",col.names=T,quote=F)
}


