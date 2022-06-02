######################################
# ARIMA for chronological genes/AS events
# usage:
# bash arima_exp.sh
# bash arima_splicing.sh

# Rscript Time_series_analysis_for_chronological.R -t $tissue -p 4 -d exp/splicing -l 0.75 > $log
######################################

library("optparse")

option_list = list(
  make_option(c("-t", "--tissue"), type="character", default=NULL, 
              help="input tissue id", metavar="character"),
  make_option(c("-p", "--threads"), type="numeric", default=10, 
              help="number of cores [default= %default]", metavar="numeric"),
  make_option(c("-d", "--dataType"), type="character", default=NULL, 
              help="number of cores [default= %default]", metavar="character"),
  make_option(c("-l", "--arima loess span"), type="numeric", default=0.75, 
              help="ARIMA loess span [default= %default]", metavar="numeric")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)


######################################
# loading packages: 
######################################
require(zoo)
require(forecast)
require(reshape2)
require(parallel)
require(dplyr)
library(data.table)

# parameters
arima.use.boxcox <- F
evalage <- NULL
arima.use.portmanteau <- T
n.cores <- opt$threads
dataType <- opt$dataType
arima.loess.span=opt$`arima loess span`
tissue_dir <- opt$tissue
tissue_dir <- gsub("\\\\","",tissue_dir)

# directory
pheno_dir <- ""
arima_dir <- paste("",dataType,"",sep = "")
TPM_dir <- ""
fitted_nonzero_dir <- paste("",dataType,"/fitted_nonzero",sep = "")
imp_PSI_dir <- ""

system(paste("mkdir ",arima_dir,sep = ""))
system(paste("mkdir ",fitted_nonzero_dir,sep = ""))

arima.selection <- function(tissue=NULL,
                            select_sex=NULL,
                            evalage=NULL,
                            use.boxcox=FALSE,
                            use.portmanteau=TRUE,
                            loess.span=0.75,
                            dataType=c("exp","splicing"),
                            n.cores=10) {
  
  select_sex_str <- ifelse(select_sex==1,"male","female")
  #
  
  ######################################
  # Data preprocessing: 
  # 1. split sex
  # 2. ordered by age: same age average
  ######################################
  cat("Building time series object\n")
  
  pheno_sex <- pheno[pheno$SEX == select_sex,]
  ages <- unique(pheno_sex$AGE)
  ages_freq <- summary(factor(pheno_sex$AGE))
  
  imp_psi_sexsel_unique <- data.frame(NA)
  for (age in ages) {
    age_freq <- ages_freq[names(ages_freq)==age]
    if (age_freq > 1){
      new <- rowMeans(imp_psi[,(colnames_imp_psi %in% pheno_sex$ID) & 
                                (colnames_imp_psi %in% pheno_sex[pheno_sex$AGE == age,]$ID)],na.rm = T)
      imp_psi_sexsel_unique <- cbind(imp_psi_sexsel_unique,new)
    }else{
      new <- imp_psi[,(colnames_imp_psi %in% pheno_sex$ID) & 
                       (colnames_imp_psi %in% pheno_sex[pheno_sex$AGE == age,]$ID)]
      imp_psi_sexsel_unique <- cbind(imp_psi_sexsel_unique,new)
    }
  }
  imp_psi_sexsel_unique <- imp_psi_sexsel_unique[,-1]
  
  colnames(imp_psi_sexsel_unique) <- paste(select_sex_str,"age",ages,sep = "_")
  
  zoo.sexsel <- zoo(t(imp_psi_sexsel_unique),order.by=ages,frequency = 1)
  ts.sexsel <- mclapply(zoo.sexsel,ts,mc.cores = n.cores)
  
  
  if (use.boxcox) {
    cat("Calculating Box-Cox Lambda using the log-likelihood method\n")
    bclambda.sexsel <- sapply(ts.sexsel,BoxCox.lambda,method="loglik",upper=10)
  } else {
    bclambda.sexsel = NULL
  }
  gc()
  
  
  ######################################
  # Builds an ARIMA model to test for trends in each time series, 
  # then selects those for which the trend is non-negligible 
  # (stationarity cannot be rejected AND no AR or MA terms are chosen)
  ######################################
  
  cat("Fitting auto-arima models on",n.cores,"CPU cores\n")
  
  arima.sexsel <- mclapply(setNames(seq_len(length(ts.sexsel)),names(ts.sexsel)),
                           function(n) auto.arima(ts.sexsel[[n]],
                                                  seasonal = F,
                                                  allowmean = T,
                                                  allowdrift = T,
                                                  lambda = bclambda.sexsel), 
                           mc.cores = n.cores)
  gc()
  
  if (use.portmanteau) {
    # portmanteau tests on RESIDUALS (Ljung-Box) - models with autocorrelated residuals are removed as well
    cat("Computing portmanteau tests on ARIMA model residuals using the Ljung-Box algorithm, on",n.cores,"CPU cores\n")
    portmanteau.sexsel <- unlist(mclapply(arima.sexsel,function(x) {
      df=sum(x$arma[1:2])
      pmin <- Box.test(x$residuals,lag=20,type="Ljung-Box",fitdf = df)$p.value
      return(pmin)
    },mc.cores = n.cores))
    cat("Proceeding to filter out",ifelse(dataType=="splicing","events","genes"),"with stationary models and autocorrelated residuals\n")
    arima.sexsel.select <- unlist(lapply(arima.sexsel,function(x) sum(x$arma[1:2])>0 & x$arma[6]>0)) & portmanteau.sexsel>0.05
  } else {
    cat("Proceeding to filter out",ifelse(dataType=="splicing","events","genes"),"with stationary models and autocorrelated residuals\n")
    arima.sexsel.select <- unlist(lapply(arima.sexsel,function(x) sum(x$arma[1:2])>0 & x$arma[6]>0))
    portmanteau.sexsel = NULL
  }
  gc()
  
  arima.sexsel.nonzero <- arima.sexsel[arima.sexsel.select]
  ts.sexsel.nonzero <- ts.sexsel[arima.sexsel.select]
  fitted.sexsel.nonzero <- t(as.data.frame(sapply(arima.sexsel.nonzero,`[[`,"fitted"),row.names = as.integer(attr(arima.sexsel.nonzero[[1]]$x,"index"))))
  cat("Out of",length(arima.sexsel),ifelse(dataType=="splicing","events,","genes"),length(ts.sexsel.nonzero),"trendy and well-fitted",
      ifelse(dataType=="splicing","splicing","genes"),"were selected for further analysis\n")
  gc()
  
  
  ######################################
  # Uses local fit (LOESS) on sex-union nonzero trends to predict 
  # chromatin values at each age in the span, for sex comparisons
  ######################################
  
  cat("Computing predicting",ifelse(dataType=="splicing","events","genes"),"using LOESS interpolation over ARIMA-fitted data on",n.cores,"cores\n")
  
  loess.sexsel.nonzero <- mclapply(arima.sexsel.nonzero,function(A) loess(as.vector(A$fitted)~attr(A$x,"index"),span = loess.span),mc.cores = n.cores)
  gc()
  predicted.sexsel.nonzero <- t(as.data.frame(mclapply(loess.sexsel.nonzero,predict,evalage,mc.cores = n.cores),row.names = evalage))
  gc()
  
  output <- list(tissue=tissue,
                 sex=select_sex,
                 ts.complete=ts.sexsel,
                 arima.complete=arima.sexsel,
                 arima.select=arima.sexsel.select,
                 loess.nonzero=loess.sexsel.nonzero,
                 fitted.nonzero=fitted.sexsel.nonzero,
                 predicted.nonzero=predicted.sexsel.nonzero,
                 BoxCox_lambda=bclambda.sexsel,
                 portmanteau.pvalues=portmanteau.sexsel)
  
  cat("Done.\n")
  return(output)
  
}

######################################
# Main script used to fit time series ARIMA models to (batch-adjusted) 
# RNA-seq gene expression and alternative splicing data
######################################

cat(tissue_dir)
cat("Starting run:\n")
cat("Fitting ARIMA model to peaks: finding trendy peaks, computing fitted and predicted data\n")

predicted.nonzero = list()
output.tmp="./"


######################################
# loading data
######################################
if(dataType == "splicing"){
  imp_psi <- read.table(paste(imp_PSI_dir,"/GTEx_v8_",tissue_dir,"_ase_psi_filter_tpm_pc_removeBE.txt",sep = ""),sep = "\t",header = T,row.names = 1)
  dim(imp_psi)
  # head(imp_psi)
}else{
  if(dataType == "exp"){
    imp_psi <- read.table(paste(TPM_dir,"/GTEx_v8_",tissue_dir,"_exp_tpm_removeBE.txt",sep = ""),sep = "\t",header = T,row.names = 1)
    imp_psi <- imp_psi[rowMeans(imp_psi)>1,]
    dim(imp_psi)
  }
}

pheno <- read.table(paste(pheno_dir,"/",tissue_dir,"_pheno.txt",sep = ""),sep = "\t")
dim(pheno) 
colnames_imp_psi <- paste(unlist(lapply(strsplit(colnames(imp_psi),".",fixed = T),function(x){x[1]})),
                          unlist(lapply(strsplit(colnames(imp_psi),".",fixed = T),function(x){x[2]})),sep = "-")
pheno <- pheno[pheno$V1 %in% colnames_imp_psi,]
rownames(pheno) <- pheno$V1
colnames(pheno) <- c("ID","SEX","AGE","RACE","BMI","DTHHRDY")
pheno <- arrange(pheno,pheno$AGE,pheno$SEX)

######################################
# start ARIMA in each sex
######################################
if (any(2 %in% pheno$SEX)) {
  cat("ARIMA fitting for females...\n============================\n")
  arima.females <- arima.selection(tissue=tissue_dir,
                                   select_sex=2,
                                   evalage=evalage,
                                   use.boxcox=arima.use.boxcox,
                                   use.portmanteau=arima.use.portmanteau,
                                   loess.span=arima.loess.span,
                                   dataType=dataType,
                                   n.cores=n.cores)
  save(arima.females,file = paste(arima_dir,"/",tissue_dir,"_","female",".RData",sep = ""))
  
  write.table(arima.females$fitted.nonzero,
              paste(fitted_nonzero_dir,"/",tissue_dir,"_","female",".txt",sep = ""),
              sep = "\t",quote = F,row.names = T,col.names = T)
  
  predicted.females.nonzero <- arima.females$predicted.nonzero
  rm("arima.females")
  gc()
  predicted.nonzero$females <- predicted.females.nonzero
}
if (any(1 %in% pheno$SEX)) {
  cat("ARIMA fitting for males...\n============================\n")
  arima.males <- arima.selection(tissue=tissue_dir,
                                 select_sex=1,
                                 evalage=evalage,
                                 use.boxcox=arima.use.boxcox,
                                 use.portmanteau=arima.use.portmanteau,
                                 loess.span=arima.loess.span,
                                 dataType=dataType,
                                 n.cores=n.cores)
  save(arima.males,file = paste(arima_dir,"/",tissue_dir,"_","male",".RData",sep = ""))
  
  write.table(arima.males$fitted.nonzero,
              paste(fitted_nonzero_dir,"/",tissue_dir,"_","male",".txt",sep = ""),
              sep = "\t",quote = F,row.names = T,col.names = T)

  predicted.males.nonzero <- arima.males$predicted.nonzero
  rm("arima.males")
  gc()
  predicted.nonzero$males <- predicted.males.nonzero
}

