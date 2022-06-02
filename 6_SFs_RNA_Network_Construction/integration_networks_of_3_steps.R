# Network construction: 
# filter overlaps
# step1: correlation between diff SFs and diff events
# step2: shRNA RBP significantly changed events overlap with diff sex
# step3: deepbind filter events with binding sites

######################################################################
# load packages:
######################################################################
library(data.table)
library(purrr)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(ggthemes)

######################################################################
# define functions:
######################################################################
fastread <- function(file,sep){
  df <- fread(file,sep = sep)
  df <- as.data.frame(df)
  rownames(df) <- df$V1
  df <- df[,-1]
  return(df)
}
bin_width <- 20
wilcox_pval_max_min_bin <- function(x){
  bin_count <- length(x)-bin_width+1
  bin_means <- c()
  for (i in 1:bin_count) {
    bin <- x[i:(i+bin_width-1)]
    bin_mean <- mean(bin)
    bin_means <- c(bin_means,bin_mean)
  }
  max_bin_pos <- which(bin_means == max(bin_means))
  min_bin_pos <- which(bin_means == min(bin_means))
  max_bin <- x[max_bin_pos:(max_bin_pos+bin_width-1)]
  min_bin <- x[min_bin_pos:(min_bin_pos+bin_width-1)]
  res <- wilcox.test(max_bin,min_bin)
  return(res$p.value)
}
ttest_pval_max_min_bin <- function(x){
  bin_count <- length(x)-bin_width+1
  bin_means <- c()
  for (i in 1:bin_count) {
    bin <- x[i:(i+bin_width-1)]
    bin_mean <- mean(bin)
    bin_means <- c(bin_means,bin_mean)
  }
  max_bin_pos <- which(bin_means == max(bin_means))
  min_bin_pos <- which(bin_means == min(bin_means))
  max_bin <- x[max_bin_pos:(max_bin_pos+bin_width-1)]
  min_bin <- x[min_bin_pos:(min_bin_pos+bin_width-1)]
  res <- t.test(max_bin,min_bin)
  return(res$p.value)
}


correlation_dir <- ""
diff_age_deepbind_dir <- ""
diff_age_pair_dir <- ""
diff_shRNA_dir <- ""
shRNA_SF_dir <- ""
tissue_dir <- ""

print(tissue_dir)
  
######################################################################
# overlap RBP: correlation & shRBP & deepbind
######################################################################
for(choose_sex in c("female","male")) {
  
  print(choose_sex)
  
  shRNA_RBP_diff_splicing_dir <- paste("/picb/rnasys2/wangsiqi/GTEx/v8/diff/SF_regulation/shRNA_SF_separately/",tissue_dir,sep = "")
  deepbind_average_dir <- paste(diff_age_deepbind_dir,"/",tissue_dir,"/deepbind_average_diff_exp_rbp_40bp_",choose_sex,sep="")
  deepbind_max_dir <- paste(diff_age_deepbind_dir,"/",tissue_dir,"/deepbind_diff_exp_rbp_40bp_",choose_sex,sep="")
  
  ######################################################################
  # step1 : filter correlation
  ######################################################################
  diff_SF_splicing_cor <- fastread(paste(correlation_dir,"/GTEx_v8_",tissue_dir,"_cor_pval_",choose_sex,".txt",sep = ""),sep = "\t")
  cor_RBP_list <- names(apply(diff_SF_splicing_cor, 1, min))[apply(diff_SF_splicing_cor, 1, min)<0.05]
  cor_events_list <- names(apply(diff_SF_splicing_cor, 2, min))[apply(diff_SF_splicing_cor, 2, min)<0.05]
  
  ######################################################################
  # step2 : filter shRNA
  ######################################################################
  diff_shRNA_RBP_id <- unlist(lapply(strsplit(list.files(diff_shRNA_dir),"_"),function(x){paste(x[1],x[2],sep = "_")}))
  
  diff_cor_shRBP_ids <- diff_shRNA_RBP_id[unlist(lapply(strsplit(diff_shRNA_RBP_id,"_"), function(x){x[1]})) %in% cor_RBP_list]
  
  if(file.exists(paste(diff_age_deepbind_dir,tissue_dir,"/",tissue_dir,"_deepbind-rbp-hs_",choose_sex,".ids",sep = ""))==F){next}
  rbp_ids <- read.delim(paste(diff_age_deepbind_dir,tissue_dir,"/",tissue_dir,"_deepbind-rbp-hs_",choose_sex,".ids",sep = ""),sep=" ",check.names = F,header = F)
  rbp_ids$V2 <- gsub("#","",rbp_ids$V2)
  
  diff_age_rbp <- rbp_ids[rbp_ids$V2 %in% unique(unlist(lapply(strsplit(diff_cor_shRBP_ids,"_"),function(x){x[1]}))),]
  if(nrow(diff_age_rbp) == 0){next;}
  
  system(paste("mkdir ",diff_age_deepbind_dir,gsub("\\)","\\\\)",gsub("\\(","\\\\(",tissue_dir)),"/deepbind_average_diff_exp_rbp_filter_plots_40bp_",choose_sex,sep = ""))
  paired_rbp_ids <- c()
  paired_rbps <- c()
  paired_events <- c()
  for (diff_cor_shRBP in unique(diff_age_rbp$V2)) {
    
    print(diff_cor_shRBP)
    diff_cor_shRBP_id <- diff_cor_shRBP_ids[which(unlist(lapply(strsplit(diff_cor_shRBP_ids,"_"),function(x){x[1]})) == diff_cor_shRBP)]
    cor_events_list_rbp <- colnames(diff_SF_splicing_cor)[as.numeric(diff_SF_splicing_cor[diff_cor_shRBP,])< 0.05]
    
    sig_diff_shRBP_events <- c()
    for (id in diff_cor_shRBP_id) {
      diff_shRBP <- read.table(paste(diff_shRNA_dir,"/",id,"_diff_splice_inter_padj.txt",sep = ""))
      sig_diff_shRBP <- diff_shRBP[diff_shRBP$pval_shRNA<0.05 & abs(diff_shRBP$deltapsi_PSI) > 0.05,]
      sig_diff_shRBP_events <- union(sig_diff_shRBP_events,sig_diff_shRBP[rownames(sig_diff_shRBP) %in% cor_events_list_rbp,]$full_id)
    }
    
    shRBP_sig_events_SE <-  sig_diff_shRBP_events[unlist(lapply(strsplit(sig_diff_shRBP_events,"@"),function(x){x[1]}))=="SE"]
    shRBP_sig_events_RI <-  sig_diff_shRBP_events[unlist(lapply(strsplit(sig_diff_shRBP_events,"@"),function(x){x[1]}))=="RI"]
    shRBP_sig_events_A5SS <-  sig_diff_shRBP_events[unlist(lapply(strsplit(sig_diff_shRBP_events,"@"),function(x){x[1]}))=="A5SS"]
    shRBP_sig_events_A3SS <-  sig_diff_shRBP_events[unlist(lapply(strsplit(sig_diff_shRBP_events,"@"),function(x){x[1]}))=="A3SS"]
    
    ######################################################################
    # step3 : filter deepbind
    ######################################################################
    if(length(shRBP_sig_events_SE)>0){
      for (cor_shRBP_GTEx_SE_event in shRBP_sig_events_SE) {
        event_symbol <- unlist(lapply(strsplit(cor_shRBP_GTEx_SE_event,"@"),function(x){x[2]}))
        
        db_res_up <- read.table(paste(deepbind_average_dir,"/SE_upstream_region/SE_upstream_region-",cor_shRBP_GTEx_SE_event,"-deepbind.txt",sep = ""),sep = "\t",header = T)
        db_res_down <- read.table(paste(deepbind_average_dir,"/SE_downstream_region/SE_downstream_region-",cor_shRBP_GTEx_SE_event,"-deepbind.txt",sep = ""),sep = "\t",header = T)
        db_res_max_up <- read.table(paste(deepbind_max_dir,"/SE_upstream_region/SE_upstream_region-",cor_shRBP_GTEx_SE_event,"-deepbind.txt",sep = ""),sep = "\t",header = T)
        db_res_max_down <- read.table(paste(deepbind_max_dir,"/SE_downstream_region/SE_downstream_region-",cor_shRBP_GTEx_SE_event,"-deepbind.txt",sep = ""),sep = "\t",header = T)
        
        rbp_db_ids <- diff_age_rbp[diff_age_rbp$V2==diff_cor_shRBP,]$V1
        for (rbp_db_id in rbp_db_ids) {
          
          if((ttest_pval_max_min_bin(db_res_max_up[,rbp_db_id])>0.05) & (ttest_pval_max_min_bin(db_res_max_down[,rbp_db_id])>0.05)){next;}
          if((max(db_res_max_up[,rbp_db_id])<1) & (max(db_res_max_down[,rbp_db_id])<1)){next;}
          print(cor_shRBP_GTEx_SE_event)
          df_plot <- data.frame(Binding_score = c(db_res_max_up[,rbp_db_id],
                                                  rep(NA,nrow(db_res_max_up)/3),
                                                  db_res_max_down[,rbp_db_id]))
          
          p <- ggplot(df_plot,aes(x=1:nrow(df_plot),y=Binding_score))+
            geom_rect(aes(xmin=-Inf, xmax=nrow(db_res_max_up), ymin=-Inf, ymax=Inf),fill='#F9F1F0',alpha = .1)+
            geom_rect(aes(xmin=nrow(df_plot)-nrow(db_res_max_down), xmax=Inf, ymin=-Inf, ymax=Inf),fill='#EEF0F5',alpha = .1)+
            geom_line()+
            geom_vline(xintercept = c(nrow(db_res_max_up),nrow(df_plot)-nrow(db_res_max_down)),linetype="dashed",color="grey20")+
            geom_vline(xintercept = c(300,nrow(df_plot)-300),linetype="dotted",color="grey20")+
            theme_few()+
            xlab("")+
            scale_x_continuous(breaks = c(1,300,nrow(db_res_max_up),nrow(df_plot)-nrow(db_res_max_down),nrow(df_plot)-300,nrow(df_plot)),
                               labels = c("-300bp","5'exon","+300bp","-300bp","3'exon","+300bp"))+
            ggtitle(paste(diff_cor_shRBP," : ",cor_shRBP_GTEx_SE_event," - ",event_symbol,sep = ""))+
            ylab("Binding Score")+
            theme(axis.text.x = element_text(size=8.5,angle=0,colour = "black",hjust = 0.5,vjust = 0.5),
                  axis.ticks.length = unit(0.05, "cm"),
                  axis.text.y = element_text(size=8.5,angle=0,colour = "black",hjust = 1,vjust = 0.5),
                  plot.title = element_text(lineheight=15, face="bold",size = 11,hjust = 0.5),
                  axis.title.y = element_text(size=9,face="bold",color="black"),
                  axis.title.x = element_text(size=0,face="bold",color="black"),
                  legend.text = element_text(size=9,color="black"))
          p
          ggsave(paste(diff_age_deepbind_dir,tissue_dir,"/deepbind_average_diff_exp_rbp_filter_plots_40bp_",choose_sex,
                       "/","RI_upstream_downstream_region","-",diff_cor_shRBP,"-",rbp_db_id,"-",cor_shRBP_GTEx_SE_event,".png",sep = ""),
                 p,width = 4.5,height = 1.8)
          paired_rbp_ids <- c(paired_rbp_ids,rbp_db_id)
          paired_rbps <- c(paired_rbps,diff_cor_shRBP)
          paired_events <- c(paired_events,cor_shRBP_GTEx_SE_event)
        }
      }
    }
    
    #### RI_upstream_region, RI_downstream_region
    if(length(shRBP_sig_events_RI)>0){
      for (cor_shRBP_GTEx_RI_event in shRBP_sig_events_RI) {
        event_symbol <- unlist(lapply(strsplit(cor_shRBP_GTEx_RI_event,"@"),function(x){x[2]}))
        
        db_res_up <- read.table(paste(deepbind_average_dir,"/RI_upstream_region/RI_upstream_region-",cor_shRBP_GTEx_RI_event,"-deepbind.txt",sep = ""),sep = "\t",header = T)
        db_res_down <- read.table(paste(deepbind_average_dir,"/RI_downstream_region/RI_downstream_region-",cor_shRBP_GTEx_RI_event,"-deepbind.txt",sep = ""),sep = "\t",header = T)
        db_res_max_up <- read.table(paste(deepbind_max_dir,"/RI_upstream_region/RI_upstream_region-",cor_shRBP_GTEx_RI_event,"-deepbind.txt",sep = ""),sep = "\t",header = T)
        db_res_max_down <- read.table(paste(deepbind_max_dir,"/RI_downstream_region/RI_downstream_region-",cor_shRBP_GTEx_RI_event,"-deepbind.txt",sep = ""),sep = "\t",header = T)
        
        rbp_db_ids <- diff_age_rbp[diff_age_rbp$V2==diff_cor_shRBP,]$V1
        for (rbp_db_id in rbp_db_ids) {
          
          if((ttest_pval_max_min_bin(db_res_max_up[,rbp_db_id])>0.05) & (ttest_pval_max_min_bin(db_res_max_down[,rbp_db_id])>0.05)){next;}
          if((max(db_res_max_up[,rbp_db_id])<1) & (max(db_res_max_down[,rbp_db_id])<1)){next;}
          print(cor_shRBP_GTEx_RI_event)
          df_plot <- data.frame(Binding_score = c(db_res_max_up[,rbp_db_id],
                                                  rep(NA,nrow(db_res_max_up)/3),
                                                  db_res_max_down[,rbp_db_id]))
          p <- ggplot(df_plot,aes(x=1:nrow(df_plot),y=Binding_score))+
            geom_rect(aes(xmin=-Inf, xmax=nrow(db_res_max_up), ymin=-Inf, ymax=Inf),fill='#F9F1F0',alpha = .1)+
            geom_rect(aes(xmin=nrow(df_plot)-nrow(db_res_max_down), xmax=Inf, ymin=-Inf, ymax=Inf),fill='#EEF0F5',alpha = .1)+
            geom_line()+
            geom_vline(xintercept = c(nrow(db_res_max_up),nrow(df_plot)-nrow(db_res_max_down)),linetype="dashed",color="grey20")+
            geom_vline(xintercept = c(300,nrow(df_plot)-300),linetype="dotted",color="grey20")+
            theme_few()+
            xlab("")+
            scale_x_continuous(breaks = c(1,300,nrow(db_res_max_up),nrow(df_plot)-nrow(db_res_max_down),nrow(df_plot)-300,nrow(df_plot)),labels = c("-300bp","5'ss","+300bp","-300bp","3'ss","+300bp"))+
            ggtitle(paste(diff_cor_shRBP," : ",cor_shRBP_GTEx_RI_event," - ",event_symbol,sep = ""))+
            ylab("Binding Score")+
            theme(axis.text.x = element_text(size=8.5,angle=0,colour = "black",hjust = 0.5,vjust = 0.5),
                  axis.ticks.length = unit(0.05, "cm"),
                  axis.text.y = element_text(size=8.5,angle=0,colour = "black",hjust = 1,vjust = 0.5),
                  plot.title = element_text(lineheight=15, face="bold",size = 11,hjust = 0.5),
                  axis.title.y = element_text(size=9,face="bold",color="black"),
                  axis.title.x = element_text(size=0,face="bold",color="black"),
                  legend.text = element_text(size=9,color="black"))
          p
          ggsave(paste(diff_age_deepbind_dir,tissue_dir,"/deepbind_average_diff_exp_rbp_filter_plots_40bp_",choose_sex,
                       "/","RI_upstream_downstream_region","-",diff_cor_shRBP,"-",rbp_db_id,"-",cor_shRBP_GTEx_RI_event,".png",sep = ""),
                 p,width = 4.5,height = 1.8)
          paired_rbp_ids <- c(paired_rbp_ids,rbp_db_id)
          paired_rbps <- c(paired_rbps,diff_cor_shRBP)
          paired_events <- c(paired_events,cor_shRBP_GTEx_RI_event)
        }
      }
    }
    
    ##### A5SS_full_region
    if(length(shRBP_sig_events_A5SS)>0){
      for (cor_shRBP_GTEx_A5SS_event in shRBP_sig_events_A5SS) {
        event_symbol <- unlist(lapply(strsplit(cor_shRBP_GTEx_A5SS_event,"@"),function(x){x[2]}))
        
        db_res <- read.table(paste(deepbind_average_dir,"/A5SS_full_region/A5SS_full_region-",cor_shRBP_GTEx_A5SS_event,"-deepbind.txt",sep = ""),sep = "\t",header = T)
        db_res_max <- read.table(paste(deepbind_max_dir,"/A5SS_full_region/A5SS_full_region-",cor_shRBP_GTEx_A5SS_event,"-deepbind.txt",sep = ""),sep = "\t",header = T)
        
        rbp_db_ids <- diff_age_rbp[diff_age_rbp$V2==diff_cor_shRBP,]$V1
        for (rbp_db_id in rbp_db_ids) {
          
          if((ttest_pval_max_min_bin(db_res_max[,rbp_db_id])>0.05) || max(db_res_max[,rbp_db_id])<1){next;}
          print(cor_shRBP_GTEx_A5SS_event)
          df_plot <- data.frame(Binding_score = db_res_max[,rbp_db_id])
          
          p <- ggplot(df_plot,aes(x=1:nrow(df_plot),y=Binding_score))+
            geom_rect(aes(xmin=-Inf, xmax=300, ymin=-Inf, ymax=Inf),fill='#F9F1F0',alpha = .1)+
            geom_rect(aes(xmin=nrow(db_res_max)-300, xmax=Inf, ymin=-Inf, ymax=Inf),fill='#EEF0F5',alpha = .1)+
            geom_line()+
            geom_vline(xintercept = c(300,nrow(db_res_max)-300),linetype="dotted",color="grey20")+
            theme_few()+
            xlab("")+
            scale_x_continuous(breaks = c(1,300,nrow(db_res_max)-300,nrow(df_plot)),labels = c("-300bp","5'ss","3'ss","+300bp"))+
            ggtitle(paste(diff_cor_shRBP," : ",cor_shRBP_GTEx_A5SS_event," - ",event_symbol,sep = ""))+
            ylab("Average Binding Score")+
            theme(axis.text.x = element_text(size=8.5,angle=0,colour = "black",hjust = 0.5,vjust = 0.5),
                  axis.ticks.length = unit(0.05, "cm"),
                  axis.text.y = element_text(size=8.5,angle=0,colour = "black",hjust = 1,vjust = 0.5),
                  plot.title = element_text(lineheight=15, face="bold",size = 11,hjust = 0.5),
                  axis.title.y = element_text(size=9,face="bold",color="black"),
                  axis.title.x = element_text(size=0,face="bold",color="black"),
                  legend.text = element_text(size=9,color="black"))
          p
          ggsave(paste(diff_age_deepbind_dir,tissue_dir,"/deepbind_average_diff_exp_rbp_filter_plots_40bp_",choose_sex,"/","A5SS_full_region","-",
                       diff_cor_shRBP,"-",rbp_db_id,"-",
                       cor_shRBP_GTEx_A5SS_event,".png",sep = ""),p,width = 3.5,height = 1.8)
          paired_rbp_ids <- c(paired_rbp_ids,rbp_db_id)
          paired_rbps <- c(paired_rbps,diff_cor_shRBP)
          paired_events <- c(paired_events,cor_shRBP_GTEx_A5SS_event)
        }
      }
    }
    
    
    ##### A3SS_full_region
    if(length(shRBP_sig_events_A3SS)>0){
      for (cor_shRBP_GTEx_A3SS_event in shRBP_sig_events_A3SS) {
        event_symbol <- unlist(lapply(strsplit(cor_shRBP_GTEx_A3SS_event,"@"),function(x){x[2]}))
        db_res <- read.table(paste(deepbind_average_dir,"/A3SS_full_region/A3SS_full_region-",cor_shRBP_GTEx_A3SS_event,"-deepbind.txt",sep = ""),sep = "\t",header = T)
        db_res_max <- read.table(paste(deepbind_max_dir,"/A3SS_full_region/A3SS_full_region-",cor_shRBP_GTEx_A3SS_event,"-deepbind.txt",sep = ""),sep = "\t",header = T)
        for (rbp_db_id in rbp_db_ids) {
          
          if((ttest_pval_max_min_bin(db_res_max[,rbp_db_id])>0.05) || max(db_res_max[,rbp_db_id])<1){next;}
          print(cor_shRBP_GTEx_A3SS_event)
          df_plot <- data.frame(Binding_score = db_res_max[,rbp_db_id])
          
          p <- ggplot(df_plot,aes(x=1:nrow(df_plot),y=Binding_score))+
            geom_rect(aes(xmin=-Inf, xmax=300, ymin=-Inf, ymax=Inf),fill='#F9F1F0',alpha = .1)+
            geom_rect(aes(xmin=nrow(db_res_max)-300, xmax=Inf, ymin=-Inf, ymax=Inf),fill='#EEF0F5',alpha = .1)+
            geom_line()+
            geom_vline(xintercept = c(300,nrow(db_res_max)-300),linetype="dotted",color="grey20")+
            theme_few()+
            xlab("")+
            scale_x_continuous(breaks = c(1,300,nrow(db_res_max)-300,nrow(df_plot)),labels = c("-300bp","5'ss","3'ss","+300bp"))+
            ggtitle(paste(diff_cor_shRBP," : ",cor_shRBP_GTEx_A3SS_event," - ",event_symbol,sep = ""))+
            ylab("Average Binding Score")+
            theme(axis.text.x = element_text(size=8.5,angle=0,colour = "black",hjust = 0.5,vjust = 0.5),
                  axis.ticks.length = unit(0.05, "cm"),
                  axis.text.y = element_text(size=8.5,angle=0,colour = "black",hjust = 1,vjust = 0.5),
                  plot.title = element_text(lineheight=15, face="bold",size = 11,hjust = 0.5),
                  axis.title.y = element_text(size=9,face="bold",color="black"),
                  axis.title.x = element_text(size=0,face="bold",color="black"),
                  legend.text = element_text(size=9,color="black"))
          p
          ggsave(paste(diff_age_deepbind_dir,tissue_dir,"/deepbind_average_diff_exp_rbp_filter_plots_40bp_",choose_sex,"/","A3SS_full_region","-",
                       diff_cor_shRBP,"-",rbp_db_id,"-",cor_shRBP_GTEx_A3SS_event,".png",sep = ""),p,
                 width = 3.5,height = 1.8)
          paired_rbp_ids <- c(paired_rbp_ids,rbp_db_id)
          paired_rbps <- c(paired_rbps,diff_cor_shRBP)
          paired_events <- c(paired_events,cor_shRBP_GTEx_A3SS_event)
        }
      }
    }
  }
  
  if(choose_sex=="male"){
    paired_male <- data.frame(paired_rbp_ids=paired_rbp_ids,
                              paired_rbps=paired_rbps,
                              paired_events=paired_events)
  }else{
    paired_female <- data.frame(paired_rbp_ids=paired_rbp_ids,
                                paired_rbps=paired_rbps,
                                paired_events=paired_events)
  }
}

######################################################################
# network filtering and integration
######################################################################

# RBP and splicing events pairs:  diff exp and AS events 
diff_sex_splice_dir <- ""
diff_sex_exp_dir <- ""
diff_splice_dir <- ""
diff_exp_dir <- ""

# filter sBASEs and sBASFs (sex-biased age-associated splicing factors)
diff_exp_results <- read.table(paste(diff_exp_dir,"/GTEx_v8_",tissue_dir,"_diff_exp_inter_MN_padj.txt",sep=""),sep = "\t",row.names = 1,header = T)
diff_sex_exp_results <- read.table(paste(diff_sex_exp_dir,"/GTEx_v8_",tissue_dir,"_diff_exp_age_MN_padj_sex_stratified.txt",sep=""),sep = "\t",row.names = 1,header = T)
diff_sex_splicing_results <- read.table(paste(diff_sex_splice_dir,"/GTEx_v8_",tissue_dir,"_diff_splice_inter_padj_sex_stratified.txt",sep=""),sep = "\t",row.names = 1,header = T)
diff_sex_splicing_results$labels <- paste(unlist(lapply(strsplit(diff_sex_splicing_results$full_id_male,"@"),function(x){paste(x[1],x[2],sep = "@")})),
                                          1:nrow(diff_sex_splicing_results),sep="-")
male_diff_exp <- diff_sex_exp_results[(diff_sex_exp_results$pval_age_class3.Old_male < 0.05) & 
                                        (abs(log2(diff_sex_exp_results$fc_TPM_age_old_male)) > log2(1.5)), ]
female_diff_exp <- diff_sex_exp_results[diff_sex_exp_results$pval_age_class3.Old_female < 0.05 &
                                          abs(log2(diff_sex_exp_results$fc_TPM_age_old_female)) > log2(1.5), ]
female_BAGES <- setdiff(female_diff_exp$symbol_female,male_diff_exp$symbol_female)
male_BAGES <- setdiff(male_diff_exp$symbol_female,female_diff_exp$symbol_female)
sUAGEs <- intersect(male_diff_exp$symbol_female,female_diff_exp$symbol_female)

# males
diff_exp_results_male_SF <- diff_sex_exp_results[diff_sex_exp_results$symbol_male %in% paired_male$paired_rbps, ]
rownames(diff_exp_results_male_SF) <- diff_exp_results_male_SF$symbol_female
diff_splicing_results_events <- diff_sex_splicing_results[diff_sex_splicing_results$full_id_male %in% paired_male$paired_events, ]
rownames(diff_splicing_results_events) <- diff_splicing_results_events$full_id_male
paired_male$diff_RBPs_log10_pval <- -log10(diff_exp_results_male_SF[paired_male$paired_rbps,]$pval_age_class3.Old_male)
paired_male$diff_RBPs_coef <- diff_exp_results_male_SF[paired_male$paired_rbps,]$coef_age_class3.Old_male
paired_male$diff_events_log10_pval <- -log10(diff_splicing_results_events[paired_male$paired_events,]$pval_age_class3.Old_male)
paired_male$labels <- diff_splicing_results_events[paired_male$paired_events,]$labels

# females
diff_exp_results_female_SF <- diff_sex_exp_results[diff_sex_exp_results$symbol_female %in% paired_female$paired_rbps, ]
rownames(diff_exp_results_female_SF) <- diff_exp_results_female_SF$symbol_female
diff_splicing_results_events <- diff_sex_splicing_results[diff_sex_splicing_results$full_id_female %in% paired_female$paired_events, ]
rownames(diff_splicing_results_events) <- diff_splicing_results_events$full_id_female
paired_female$diff_RBPs_log10_pval <- -log10(diff_exp_results_female_SF[paired_female$paired_rbps,]$pval_age_class3.Old_female)
paired_female$diff_RBPs_coef <- diff_exp_results_female_SF[paired_female$paired_rbps,]$coef_age_class3.Old_female
paired_female$diff_events_log10_pval <- -log10(diff_splicing_results_events[paired_female$paired_events,]$pval_age_class3.Old_female)
paired_female$labels <- diff_splicing_results_events[paired_female$paired_events,]$labels
dim(paired_male)
dim(paired_female)

# labeled weight of edges
diff_SF_splicing_cor <- fastread(paste(correlation_dir,"/GTEx_v8_",tissue_dir,"_cor_pval_","female",".txt",sep = ""),sep = "\t")
edge_weights <- c()
for (i in 1:nrow(paired_female)) {
  edge_weight <- -log10(diff_SF_splicing_cor[paired_female[i,]$paired_rbps,paired_female[i,]$paired_events])
  edge_weights <- c(edge_weights,edge_weight)
}
paired_female$edges_weight <- edge_weights
diff_SF_splicing_cor <- fastread(paste(correlation_dir,"/GTEx_v8_",tissue_dir,"_cor_pval_","male",".txt",sep = ""),sep = "\t")
edge_weights <- c()
for (i in 1:nrow(paired_male)) {
  edge_weight <- -log10(diff_SF_splicing_cor[paired_male[i,]$paired_rbps,paired_male[i,]$paired_events])
  edge_weights <- c(edge_weights,edge_weight)
}
paired_male$edges_weight <- edge_weights
diff_paired_male <- paired_male[paired_male$paired_rbps %in% male_diff_exp$symbol_male,]
diff_paired_female <- paired_female[paired_female$paired_rbps %in% female_diff_exp$symbol_female,]
write.table(diff_paired_male,
            paste(diff_age_pair_dir,"/",tissue_dir,"_age_diff_RBP_events_pair_male.txt",sep=""),
            quote = F,sep = "\t",row.names = F)
write.table(diff_paired_female,
            paste(diff_age_pair_dir,"/",tissue_dir,"_age_diff_RBP_events_pair_female.txt",sep=""),
            quote = F,sep = "\t",row.names = F)

# filtering redundant regulations
paired_male_unique <-  diff_paired_male[,-1] %>% unique()
paired_female_unique <-  diff_paired_female[,-1] %>% unique()
write.table(paired_male_unique,
            paste(diff_age_pair_dir,"/",tissue_dir,"_age_diff_RBP_events_pair_male_unique.txt",sep=""),
            quote = F,sep = "\t",row.names = F)
write.table(paired_female_unique,
            paste(diff_age_pair_dir,"/",tissue_dir,"_age_diff_RBP_events_pair_female_unique.txt",sep=""),
            quote = F,sep = "\t",row.names = F)

# regulations of sBASEs
male_diff_splicing <- diff_sex_splicing_results[(diff_sex_splicing_results$pval_age_class3.Old_male < 0.05) & 
                                                  (abs(diff_sex_splicing_results$deltapsi_PSI_age_old_male) > 0.05), ]
female_diff_splicing <- diff_sex_splicing_results[(diff_sex_splicing_results$pval_age_class3.Old_female < 0.05) & 
                                                    (abs(diff_sex_splicing_results$deltapsi_PSI_age_old_female) > 0.05), ]
female_BASES <- setdiff(female_diff_splicing$full_id_female,male_diff_splicing$full_id_male)
male_BASES <- setdiff(male_diff_splicing$full_id_female,female_diff_splicing$full_id_male)
paired_male_unique_sBASEs <- paired_male_unique[paired_male_unique$paired_events %in% male_BASES,]
paired_female_unique_sBASEs <- paired_female_unique[paired_female_unique$paired_events %in% female_BASES,]
paired_male_unique_sBASEs$edge <- rep("male",nrow(paired_male_unique_sBASEs))
paired_female_unique_sBASEs$edge <- rep("female",nrow(paired_female_unique_sBASEs))
paired_male_unique_sBASEs$rbp_class <- ifelse(paired_male_unique_sBASEs$paired_rbps %in% sUAGEs,"sUAGEs",
                                              ifelse(paired_male_unique_sBASEs$paired_rbps %in% female_BAGES,"female_BAGES",
                                                     ifelse(paired_male_unique_sBASEs$paired_rbps %in% male_BAGES,"male_BAGES","None")))
paired_female_unique_sBASEs$rbp_class <- ifelse(paired_female_unique_sBASEs$paired_rbps %in% sUAGEs,"sUAGEs",
                                                ifelse(paired_female_unique_sBASEs$paired_rbps %in% female_BAGES,"female_BAGES",
                                                       ifelse(paired_female_unique_sBASEs$paired_rbps %in% male_BAGES,"male_BAGES","None")))
write.table(paired_male_unique_sBASEs,
            paste(diff_age_pair_dir,"/",tissue_dir,"_age_diff_RBP_events_pair_male_sBASEs_unique.txt",sep=""),
            quote = F,sep = "\t",row.names = F)
write.table(paired_female_unique_sBASEs,
            paste(diff_age_pair_dir,"/",tissue_dir,"_age_diff_RBP_events_pair_female_sBASEs_unique.txt",sep=""),
            quote = F,sep = "\t",row.names = F)

# sex-specific regulations
paired_male_unique_sBASEs_sp <- paired_male_unique_sBASEs[paired_male_unique_sBASEs$paired_rbps %in% setdiff(paired_male_unique_sBASEs$paired_rbps,paired_female_unique_sBASEs$paired_rbps),]
paired_female_unique_sBASEs_sp <- paired_female_unique_sBASEs[paired_female_unique_sBASEs$paired_rbps %in% setdiff(paired_female_unique_sBASEs$paired_rbps,paired_male_unique_sBASEs$paired_rbps),]
write.table(paired_male_unique_sBASEs_sp,
            paste(diff_age_pair_dir,"/",tissue_dir,"_age_diff_RBP_events_pair_male_sBASEs_sp_unique.txt",sep=""),
            quote = F,sep = "\t",row.names = F)
write.table(paired_female_unique_sBASEs_sp,
            paste(diff_age_pair_dir,"/",tissue_dir,"_age_diff_RBP_events_pair_female_sBASEs_sp_unique.txt",sep=""),
            quote = F,sep = "\t",row.names = F)

# merge networks in females and males
paired_male_unique_sBASEs$edge <- rep("male",nrow(paired_male_unique_sBASEs))
paired_female_unique_sBASEs$edge <- rep("female",nrow(paired_female_unique_sBASEs))
paired_merge_unique_sBASEs <- rbind(paired_male_unique_sBASEs,paired_female_unique_sBASEs)
paired_merge_unique_sBASEs$rbp_class <- ifelse(paired_merge_unique_sBASEs$paired_rbps %in% sUAGEs,"sUAGEs",
                                               ifelse(paired_merge_unique_sBASEs$paired_rbps %in% female_BAGES,"female_BAGES",
                                                      ifelse(paired_merge_unique_sBASEs$paired_rbps %in% male_BAGES,"male_BAGES","None")))
write.table(paired_merge_unique_sBASEs,
            paste(diff_age_pair_dir,"/",tissue_dir,"_age_diff_RBP_events_pair_merge_sBASEs_unique.txt",sep=""),
            quote = F,sep = "\t",row.names = F)

