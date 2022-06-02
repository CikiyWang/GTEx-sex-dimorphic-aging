######################################
# breakpoint analysis for subsampled genes for detecting ADGs
######################################

setwd("/picb/rnasys2/wangsiqi/GTEx/v8/time_series/breakpoint")

######################################
# loading packages: 
######################################
require(reshape2)
require(ggplot2)
require(tibble)
require(dplyr)
require(mclust)
require(MKmisc)
library("RColorBrewer")
library("ggsci")
require(parallel)
library("ggthemes")


######################################
# set parameters: 
######################################
winspan = 5:15 # Range of window sizes tested in moving-window analyses
mininterval = 5 # minimum number of years separating two breakpoints. Breakpoints occuring within this interval are averaged by age
p.thresh = 0.05 # P-value cutoff to regard a local minumum P-value as significant
range.thresh = 0.1 # threshhold on the proportion of maximum sum of  pvalue that is spanned between a minimum and the nearest maximum. Minima separated from maxima by quantities larger than this value are considered breakpoints
stepsize = 1 # step size used in moving window analysis
bwspan = seq(0.25,0.75,0.05)

dataType <- "exp"

n.cores <- 10

fitted_nonzero_dir <- paste("",dataType,"/merge_sUAGEs_pc/sUAGEs_ts_matrix",sep = "")
tissue_dirs_male <- gsub("_male_sUAGEs_ts.txt","",list.files(fitted_nonzero_dir,"_male_sUAGEs_ts.txt"))
tissue_dirs_female <- gsub("_female_sUAGEs_ts.txt","",list.files(fitted_nonzero_dir,"_female_sUAGEs_ts.txt"))
tissue_dirs <- intersect(tissue_dirs_male,tissue_dirs_female)

for (tissue_dir in tissue_dirs) {
  print(tissue_dir)
  #tissue_dir <- "Adipose-Subcutaneous"
  
  pc_ENSG_symbol <- read.table("gencode.v36.protein_coding_gene_ENSG_symbol.bed",sep = "\t")
  
  male <- read.table(paste(fitted_nonzero_dir,"/",tissue_dir,"_male_sUAGEs_ts.txt",sep = ""),
                     sep = "\t",header = T,row.names = 1)
  colnames(male) <- gsub("X","",colnames(male))
  
  female <- read.table(paste(fitted_nonzero_dir,"/",tissue_dir,"_female_sUAGEs_ts.txt",sep = ""),
                       sep = "\t",header = T,row.names = 1)
  colnames(female) <- gsub("X","",colnames(female))
  
  male <- male[rownames(male) %in% pc_ENSG_symbol$V5,]
  female <- female[rownames(female) %in% pc_ENSG_symbol$V5,]
  
  if(nrow(male)<10 || nrow(female)<10){next}
  
  clustered.datas_male <- setNames(list(male),
                                   "merge")
  clustered.datas_female <- setNames(list(female),
                                     "merge")
  clustered.datas <- list(males=clustered.datas_male,
                          females=clustered.datas_female)
  testage.list <- list(males=as.numeric(colnames(male)),
                       females=as.numeric(colnames(female)))
  
  age_ts_region <- function(x,min_winspan){
    # seqs <- seq(min(x)+min_winspan,max(x)-min_winspan,by=1)
    # seqs <- seq(min(x),max(x),by=1)
    ranges <- c()
    for (i in x[which(x >= min(x)+min_winspan & x <= max(x)-min_winspan)]) {
      len_left <- length(intersect(x,(i - min_winspan):(i)))
      len_right <- length(intersect(x,(i):(i + min_winspan)))
      ranges <- c(ranges,min(len_left,len_right))
    }
    return(min(ranges))
  }
  
  print(age_ts_region(testage.list$males,4))
  print(age_ts_region(testage.list$females,4))
  
  if((age_ts_region(testage.list$males,4)<2) || (age_ts_region(testage.list$females,4)<2)){
    cat("Sample sizes not enough; Programming breaks...");next}
  
  bk_dir <- paste("",dataType,"/merge_sUAGEs_pc/",tissue_dir,"_merge_sUAGEs_cut_0.2",sep = "")
  tissue_dir_mkdir <- gsub("\\)","\\\\)",gsub("\\(","\\\\(",tissue_dir))
  mk_bk_dir <- paste("",dataType,"/merge_sUAGEs_pc/",tissue_dir_mkdir,"_merge_sUAGEs_cut_0.2",sep = "")
  system(paste("mkdir ",mk_bk_dir))
  
  setwd(bk_dir)
  
  
  desample_count_female <- round(0.2*nrow(female))
  desample_count_male <- round(0.2*nrow(male))
  
  df_sampling_male <- data.frame()
  df_sampling_female <- data.frame()
  set.seed(1001)
  
  l <- 200
  for (i in 1:l) {
    df_sampling_male <- rbind(df_sampling_male,sample(1:nrow(male),desample_count_male))
    df_sampling_female <- rbind(df_sampling_female,sample(1:nrow(female),desample_count_female))
  }
  rownames(df_sampling_male) <- paste("round",1:l,sep = "_")
  rownames(df_sampling_female) <- paste("round",1:l,sep = "_")
  
  write.table(df_sampling_male,
              paste(bk_dir,"/",tissue_dir,"_desampling_male_cut_0.2.txt",sep = ""),
              sep = "\t",row.names = T,col.names = T,quote = F)
  write.table(df_sampling_female,
              paste(bk_dir,"/",tissue_dir,"_desampling_female_cut_0.2.txt",sep = ""),
              sep = "\t",row.names = T,col.names = T,quote = F)
  
  for (i in 1:l) {
    print(i)
    
    clustered.datas_male <- setNames(list(male[-as.numeric(df_sampling_male[i,]),]),"merge")
    clustered.datas_female <- setNames(list(female[-as.numeric(df_sampling_female[i,]),]),"merge")
    clustered.datas <- list(males=clustered.datas_male,
                            females=clustered.datas_female)
    
    testage.list <- list(males=as.numeric(colnames(male)),
                         females=as.numeric(colnames(female)))
    
    cat("Start calculate breakpoints_bysex ...")
    breakpoints_bysex <- lapply(list(males="males",females="females"), function(testsamps) {
      print(testsamps)
      testage=testage.list[[testsamps]]
      print(testage)
      clustered.data_sex <- clustered.datas[[testsamps]]
      nclust = length(clustered.data_sex)
      breakpoints_bycluster <- lapply(setNames(as.list(winspan),winspan), function(ws) {
        # Breakpoints are identified for each window span separately
        cat("Extracting breakpoints for winspan=",ws,"\n")
        offset=c(min(testage)+ws,max(testage)-ws)
        testage.offset <- seq(offset[1],offset[2],stepsize)
        lapply(clustered.data_sex, function(Z) {
          # cat(Z$Cluster[1],"\n")
          
          Z <- t(apply(Z, 1, function(y) (y-mean(y))/sd(y)))
          
          pvalues <- do.call(rbind,mclapply(bwspan, function(bw) {
            # For a given window span, breakpoints result from the integration of smoother p-value time trends across
            # multiple smoothing bandwidths. These p-values estimate the magnitude of the difference between ages prior
            # and posterior to each age point (where window span defines the sampled age period, and therefore sample size
            # and statistical power of the test)
            cat("Bandwidth=",bw,"\n")
            x <- t(sapply(setNames(as.list(testage.offset),as.numeric(testage.offset)), function(x) {
              # print(x)
              agelo <- x-ws
              samplo <- testage[testage>=agelo & testage<=x]
              agehi <- x+ws
              samphi <- testage[testage>=x & testage<=agehi]
              Y <- t(cbind(Z[,colnames(Z) %in% samplo],Z[,colnames(Z) %in% samphi]))
              Y <- Y[, colnames(Y) %in% names(which(apply(Y, 2, sd)>0))]
              V <- prcomp(Y,scale. = T) # Dimensionality reduction prior to adjacent-window comparison
              Y$win <- c(rep("low",ncol(as.matrix(Z[,colnames(Z) %in% samplo]))),rep("high",ncol(as.matrix(Z[,colnames(Z) %in% samphi]))))
              wintest <- summary(manova(V$x[,1:3] ~ Y$win))$stats[1,6] # Input first 3 PCs to MANOVA for between-window comparison
              
              # add sample sizes to interpret p values under the light of power to detect differences
              output <- c(p=wintest,n_lo=ncol(Z[,colnames(Z) %in% samplo]),n_hi=ncol(Z[,colnames(Z) %in% samphi]))
              return(output)
            })) %>% 
              data.frame(.,row.names = NULL,check.names = F) %>% 
              mutate(age=testage.offset)
            # P-value smoothing in preparation for detection of local minima
            y <- loess(data=x,formula = p ~ age,span=bw)
            z <- data.frame(age=y$x,observed_p=y$y,fitted_p=y$fitted,n_low=x$n_lo,n_high=x$n_hi,bandwidth=bw,winspan=ws)
            return(z)
          },mc.cores = n.cores))
          
          # Truncating resulting p-values to (0,1] interval
          pvalues <- pvalues %>%
            mutate(fitted_p=ifelse(fitted_p>1,1,
                                   ifelse(fitted_p<=0,min(.$fitted_p[.$fitted_p>0]),fitted_p)))
          
          # Detection of p-value maxima/minima is done on smnoothed p-values integrated by each evaluation age point
          # At each age, p-values calculated at all bandwidths for a given window span are integrated using Fisher's 
          # method, and the resulting values are smoothed again. Local minima and maxima are located in this function
          # using numerical differentiation. A parametric "significant" value is also calculated for each age point
          # by comparing the smoothed Fisher's P value to its null Chi-squared expectation.
          
          breakpoints = pvalues %>% 
            group_by(age) %>% 
            summarize(sum_fitted_p=-2*sum(log(fitted_p)), # Fisher's method for P-value combination
                      sum_observed_p=-2*sum(log(observed_p)),
                      mean_fitted_p=mean(-log10(fitted_p)), # Fisher's method for P-value combination
                      mean_observed_p=mean(-log10(observed_p))) %>%
            ungroup() %>% 
            mutate(mean_refitted_p = loess(mean_fitted_p ~ age, data = data.frame(.),span = 0.5)$fitted) %>% # 
            mutate(sum_refitted_p = loess(sum_fitted_p ~ age, data = data.frame(.),span = 0.5)$fitted)
          return(list(pvalues=pvalues,breakpoints=breakpoints))
        })
      })
    })
    
    # Diagnostic plots: p value trajectories marking breakpoints, breakpoints with variation (error bars)
    cat("Start calculate allBreaks ...")
    allBreaks <- lapply(breakpoints_bysex,
                        function(X) {
                          nclust = length(X[[1]])
                          lapply(setNames(as.list(seq_len(nclust)),seq_len(nclust)),
                                 function(k)
                                   lapply(lapply(X,
                                                 function(S) lapply(S,`[[`,"breakpoints")),`[[`,k))
                        })
    
    # output mean_refitted_p in at age point 
    res <- do.call(rbind,lapply(names(allBreaks),
                                function(sx)
                                  do.call(rbind, lapply(names(allBreaks[[sx]]),
                                                        function(k) do.call(rbind,lapply(names(allBreaks[[sx]][[k]]),
                                                                                         function(s) allBreaks[[sx]][[k]][[s]] %>%
                                                                                           dplyr::select(age,mean_refitted_p,sum_refitted_p) %>%
                                                                                           mutate(sex=sx,cluster=k,winspan=s))))))) %>%
      group_by(sex,cluster,age) %>%
      summarize(median_win_refitted_p=median(mean_refitted_p),
                mean_win_refitted_p=mean(mean_refitted_p),
                mean_win_sum_refitted_p=mean(sum_refitted_p)) %>%
      ungroup()
    
    all_breakpoints <- data.frame()
    for (sex in c("females","males")) {
      res_sex <- res[res$sex == sex,]

      breakpoints = res_sex %>% 
        mutate(mean_win_refitted_p = loess(mean_win_refitted_p ~ age, data = data.frame(.),span = 0.25)$fitted) %>%
        mutate(mean_win_sum_refitted_p = loess(mean_win_sum_refitted_p ~ age, data = data.frame(.),span = 0.25)$fitted) %>%
        mutate(isMinimum=age %in% .[c(which(c(diff(sign(diff(.$mean_win_refitted_p))),FALSE,FALSE)==2)+1,
                                      which(c(FALSE,diff(.$mean_win_refitted_p)<0) & age==max(.$age)),
                                      which(c(diff(.$mean_win_refitted_p)>0,FALSE) & age==min(.$age))),]$age,
               isMaximum=age %in% .[c(which(c(diff(sign(diff(.$mean_win_refitted_p))),FALSE,FALSE)==-2)+1,
                                      which(c(diff(.$mean_win_refitted_p)<0,FALSE) & age==min(.$age)),
                                      which(c(FALSE,diff(.$mean_win_refitted_p)>0) & age==max(.$age))),]$age)
      if (!any(breakpoints$isMinimum)) {
        breakpoints <- breakpoints %>%
          mutate(isMinimum=(mean_win_refitted_p==min(.$mean_win_refitted_p)[1]))
      }
      if (!any(breakpoints$isMaximum)) {
        breakpoints <- breakpoints %>%
          mutate(isMaximum=(mean_win_refitted_p==max(.$mean_win_refitted_p)[1]))
      }
      
      # As a second "significance" criterion, for each minimum we find the distance to the nearest *previous (in time)* 
      # maximum: drops that exceed a pre-defined threshold (range.thresh) are considered significant as well
      MinMax <- list(Min=(breakpoints %>% filter(isMinimum))$age,Max=(breakpoints %>% filter(isMaximum))$age)
      MinMaxMx <- outer(seq(length(MinMax$Min)),
                        seq(length(MinMax$Max)),
                        Vectorize(function(i,j) abs(MinMax$Min[i]-MinMax$Max[j])))
      MinMaxMatch <- data.frame(age=MinMax$Min[apply(MinMaxMx,2,which.min)]) %>%
        inner_join(breakpoints %>% dplyr::select(age,mean_win_refitted_p),by="age") %>%
        rename_(MinAge="age",min_mean_refitted_p="mean_win_refitted_p") %>%
        cbind(data.frame(age=MinMax$Max))
      
      # Finally, a breakpoint is defined for a specific window span for age points which satisfy the two significance
      # conditions calculated above
      breakpoints <- breakpoints %>% 
        left_join(MinMaxMatch,by="age") %>% 
        mutate(diffToMin=ifelse(!is.na(MinAge),(mean_win_refitted_p-min_mean_refitted_p)/max(.$mean_win_refitted_p),0)) %>% 
        mutate(isBreak=diffToMin>=range.thresh)
      
      # Breakpoint post-processing:
      # Collapsing breakpoints that are too close (within ws years, i.e. whose difference is lower than the resolution of a test given window span)
      # Algorithm: sort breakpoints by age and work from the oldest collapsing any breakppoint within the restricting interval, 
      # keeping only the one with the lowest Pvalue
      ws <- 5
      breaks <- subset(breakpoints,isBreak)
      nBreaks = nrow(breaks)
      if (nBreaks>1) {
        mBreaks = 0
        while(nBreaks!=mBreaks) {
          mBreaks = nBreaks
          brlags <- cbind(breaks$age,lag(breaks$age))
          breaks <- subset(breaks,!age %in% brlags[last(which(abs(apply(brlags,1,diff))<ws)),]) %>%
            rbind(subset(breaks,age %in% brlags[last(which(abs(apply(brlags,1,diff))<ws)),]) %>% arrange(mean_win_refitted_p) %>% slice_(1))
          nBreaks = nrow(breaks)
        }
      } 
      
      breakpoints <- rbind(subset(breakpoints,!age %in% breaks$age) %>% mutate(isBreak=F),breaks) %>%
        arrange(age)
      
      breakpoints$isBreak[1] <- FALSE
      breakpoints$isBreak[nrow(breakpoints)] <- FALSE
      breaks <- subset(breakpoints,isBreak)
      nBreaks = nrow(breaks)
      
      if (nBreaks==0){
        max_pos <- which(breakpoints$mean_win_refitted_p == max(breakpoints$mean_win_refitted_p))
        breakpoints$isBreak[max_pos] <- TRUE
      }
      
      all_breakpoints <- rbind(all_breakpoints,breakpoints)
    }
    
    write.table(as.data.frame(all_breakpoints),
                paste(bk_dir,"/",tissue_dir,"_dis_0.1_loess_0.25_mean_all_cut_0.2_round",i,".txt",sep = ""),
                sep = "\t",row.names = T,col.names = T,quote = F)
    
    rm(breakpoints_bysex)
  }
}


