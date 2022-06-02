######################################################################
# step3 : filter diff AS events with RBP binding signals by deepbind
######################################################################

library("data.table")
library("dplyr")
require("parallel")

fastread <- function(file,sep){
  df <- fread(file,sep = sep)
  df <- as.data.frame(df)
  rownames(df) <- df$V1
  df <- df[,-1]
  return(df)
}

diff_splice_dir <- ""
diff_age_deepbind_dir <- ""
diff_shRNA_dir <- ""
shRNA_SF_dir <- ""
cor_dir <- ""
window <- 40
region <- 300
tissue_dir <- ""

SF <- read.table("GO_RNA_SPLICING_SFs.txt",sep = "\t",header = F)$V1
rbp_ids <- read.delim("deepbind-rbp-hs.ids",sep = " ",check.names = F,header = F)
rbp_ids$V2 <- gsub("#","",rbp_ids$V2)
rbp_ids <- rbind(rbp_ids,t(matrix(c("D00210.001","RBFOX2")))) # we consider the similar binding motifs in RBFOX family 

print(tissue_dir)
diff_splicing_results <- read.table(paste(diff_splice_dir,"/GTEx_v8_",tissue_dir,"_diff_splice_inter_padj_sex_stratified.txt",sep=""),sep = "\t",row.names = 1,header = T)
diff_splicing_age_female <- diff_splicing_results[(diff_splicing_results$pval_age_class3.Old_female < 0.05) & (abs(diff_splicing_results$deltapsi_PSI_age_old_female) > 0.05 ),]
diff_splicing_age_male <- diff_splicing_results[(diff_splicing_results$pval_age_class3.Old_male < 0.05) & (abs(diff_splicing_results$deltapsi_PSI_age_old_male) > 0.05 ),]

# filter annotation files with spliced sites
anno <- data.frame(full_id = diff_splicing_results$full_id_female)
rownames(anno) <- anno$full_id
anno$geneSymbol <- unlist(lapply(strsplit(anno$full_id,"@"),function(x){x[2]}))
anno$chr <- unlist(lapply(strsplit(anno$full_id,"@"),function(x){x[3]}))
anno$event_type <- unlist(lapply(strsplit(anno$full_id,"@"),function(x){x[1]}))
anno$strand <- unlist(lapply(strsplit(anno$full_id,"@"),function(x){x[4]}))
diff_splicing_anno_SE_female <- anno[anno$full_id %in% diff_splicing_age_female$full_id_female & anno$event_type == "SE",]
diff_splicing_anno_RI_female <- anno[anno$full_id %in% diff_splicing_age_female$full_id_female & anno$event_type == "RI",]
diff_splicing_anno_A5SS_female <- anno[anno$full_id %in% diff_splicing_age_female$full_id_female & anno$event_type == "A5SS",]
diff_splicing_anno_A3SS_female <- anno[anno$full_id %in% diff_splicing_age_female$full_id_female & anno$event_type == "A3SS",]
diff_splicing_anno_SE_male <- anno[anno$full_id %in% diff_splicing_age_male$full_id_male & anno$event_type == "SE",]
diff_splicing_anno_RI_male <- anno[anno$full_id %in% diff_splicing_age_male$full_id_male & anno$event_type == "RI",]
diff_splicing_anno_A5SS_male <- anno[anno$full_id %in% diff_splicing_age_male$full_id_male & anno$event_type == "A5SS",]
diff_splicing_anno_A3SS_male <- anno[anno$full_id %in% diff_splicing_age_male$full_id_male & anno$event_type == "A3SS",]


######################################################################
# make bed files in adjacent spliced regions
######################################################################
# SE: upstream intron 300bp, whole exon, downsteam intron 300bp
# RI: upstream exon 300bp, upstream intron 200 bp; downstream intron 300bp, downstream exon 300bp
# A5SS: 1st ss upsteam exon 300bp, 1st ss to 2nd ss, 2nd ss downstream intron 300bp
# A3SS: 1st ss upsteam intron 300bp, 1st ss to 2nd ss, 2nd ss downstream exon 300bp
# MXE: 1st exon upstream 300bp to downstream 300bp; 2nd exon upstream 300bp to downstream 300bp

# one event one .fa
system(paste("mkdir ",diff_age_deepbind_dir,gsub("\\)","\\\\)",gsub("\\(","\\\\(",tissue_dir)),sep = ""))
system(paste("mkdir ",diff_age_deepbind_dir,gsub("\\)","\\\\)",gsub("\\(","\\\\(",tissue_dir)),"/bed6_300bp_female",sep = ""))
system(paste("mkdir ",diff_age_deepbind_dir,gsub("\\)","\\\\)",gsub("\\(","\\\\(",tissue_dir)),"/fa_300bp_female",sep = ""))
system(paste("mkdir ",diff_age_deepbind_dir,gsub("\\)","\\\\)",gsub("\\(","\\\\(",tissue_dir)),"/seq_window_40bp_female",sep = ""))
system(paste("mkdir ",diff_age_deepbind_dir,gsub("\\)","\\\\)",gsub("\\(","\\\\(",tissue_dir)),"/bed6_300bp_male",sep = ""))
system(paste("mkdir ",diff_age_deepbind_dir,gsub("\\)","\\\\)",gsub("\\(","\\\\(",tissue_dir)),"/fa_300bp_male",sep = ""))
system(paste("mkdir ",diff_age_deepbind_dir,gsub("\\)","\\\\)",gsub("\\(","\\\\(",tissue_dir)),"/seq_window_40bp_male",sep = ""))


for (choose_sex in c("female","male")) {
  print(choose_sex)
  if(choose_sex == "female"){
    diff_splicing_anno_SE <- diff_splicing_anno_SE_female
    diff_splicing_anno_RI <- diff_splicing_anno_RI_female
    diff_splicing_anno_A5SS <- diff_splicing_anno_A5SS_female
    diff_splicing_anno_A3SS <- diff_splicing_anno_A3SS_female
  }
  if(choose_sex == "male"){
    diff_splicing_anno_SE <- diff_splicing_anno_SE_male
    diff_splicing_anno_RI <- diff_splicing_anno_RI_male
    diff_splicing_anno_A5SS <- diff_splicing_anno_A5SS_male
    diff_splicing_anno_A3SS <- diff_splicing_anno_A3SS_male
  }
  
  # filter existed RBPs in previous 2 steps first
  cors_pval <- fastread(paste(cor_dir,"/GTEx_v8_",tissue_dir,"_cor_pval_",choose_sex,".txt",sep = ""),sep = "\t")
  rbp_cor_pval <- apply(cors_pval, 1, min)
  events_cor_pval <- apply(cors_pval, 2, min)
  cor_RBP_symbols <- names(rbp_cor_pval[rbp_cor_pval<0.05])
  cor_events_ids <- names(events_cor_pval[events_cor_pval<0.05])
  
  diff_shRNA_RBP_id <- unlist(lapply(strsplit(list.files(diff_shRNA_dir),"_"),function(x){paste(x[1],x[2],sep = "_")}))
  diff_cor_shRBP_ids <- diff_shRNA_RBP_id[unlist(lapply(strsplit(diff_shRNA_RBP_id,"_"), function(x){x[1]})) %in% cor_RBP_symbols]
  diff_age_rbp <- rbp_ids[rbp_ids$V2 %in% unique(unlist(lapply(strsplit(diff_cor_shRBP_ids,"_"),function(x){x[1]}))),]
  if(nrow(diff_age_rbp) == 0){next;}
  rbp_ids_tissue <- data.frame(diff_age_rbp$V1,paste("#",diff_age_rbp$V2,sep=""))
  write.table(rbp_ids_tissue,paste(diff_age_deepbind_dir,tissue_dir,"/",tissue_dir,"_deepbind-rbp-hs_",choose_sex,".ids",sep = ""),quote = F,sep = " ",row.names = F,col.names = F)
  
  # filter existed AS events in previous 2 steps first
  shRBP_sig_events_SEs <- c()
  shRBP_sig_events_RIs <- c()
  shRBP_sig_events_A5SSs <- c()
  shRBP_sig_events_A3SSs <- c()
  for (diff_cor_shRBP in unique(diff_age_rbp$V2)) {
    diff_cor_shRBP_id <- diff_cor_shRBP_ids[which(unlist(lapply(strsplit(diff_cor_shRBP_ids,"_"),function(x){x[1]})) == diff_cor_shRBP)]
    sig_diff_shRBP_events <- c()
    for (id in diff_cor_shRBP_id) {
      diff_shRBP <- read.table(paste(diff_shRNA_dir,"/",id,"_diff_splice_inter_padj.txt",sep = ""))
      sig_diff_shRBP <- diff_shRBP[diff_shRBP$pval_shRNA<0.05 & abs(diff_shRBP$deltapsi_PSI) > 0.05,]
      sig_diff_shRBP_events <- union(sig_diff_shRBP_events,sig_diff_shRBP[rownames(sig_diff_shRBP) %in% cor_events_ids,]$full_id)
    }
    shRBP_sig_events_SE <-  sig_diff_shRBP_events[unlist(lapply(strsplit(sig_diff_shRBP_events,"@"),function(x){x[1]}))=="SE"]
    shRBP_sig_events_RI <-  sig_diff_shRBP_events[unlist(lapply(strsplit(sig_diff_shRBP_events,"@"),function(x){x[1]}))=="RI"]
    shRBP_sig_events_A5SS <-  sig_diff_shRBP_events[unlist(lapply(strsplit(sig_diff_shRBP_events,"@"),function(x){x[1]}))=="A5SS"]
    shRBP_sig_events_A3SS <-  sig_diff_shRBP_events[unlist(lapply(strsplit(sig_diff_shRBP_events,"@"),function(x){x[1]}))=="A3SS"]
    
    shRBP_sig_events_SEs <- union(shRBP_sig_events_SEs,shRBP_sig_events_SE)
    shRBP_sig_events_RIs <- union(shRBP_sig_events_RIs,shRBP_sig_events_RI)
    shRBP_sig_events_A5SSs <- union(shRBP_sig_events_A5SSs,shRBP_sig_events_A5SS)
    shRBP_sig_events_A3SSs <- union(shRBP_sig_events_A3SSs,shRBP_sig_events_A3SS)
  }
  diff_splicing_anno_SE_sel <- diff_splicing_anno_SE[diff_splicing_anno_SE$full_id %in% shRBP_sig_events_SEs,]
  diff_splicing_anno_RI_sel <- diff_splicing_anno_RI[diff_splicing_anno_RI$full_id %in% shRBP_sig_events_RIs,]
  diff_splicing_anno_A5SS_sel <- diff_splicing_anno_A5SS[diff_splicing_anno_A5SS$full_id %in% shRBP_sig_events_A5SSs,]
  diff_splicing_anno_A3SS_sel <- diff_splicing_anno_A3SS[diff_splicing_anno_A3SS$full_id %in% shRBP_sig_events_A3SSs,]
  
  if(nrow(diff_splicing_anno_SE_sel)>0){
    diff_splicing_anno_SE_sel$splice_out <- unlist(lapply(strsplit(diff_splicing_anno_SE_sel$full_id,"@"),function(x){x[6]}))
    diff_splicing_anno_SE_sel$start <-  as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(diff_splicing_anno_SE_sel$splice_out,"^",fixed = T),
                                                                                        function(x){x[2]})),"-"),function(x){x[1]})))
    diff_splicing_anno_SE_sel$end <-  as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(diff_splicing_anno_SE_sel$splice_out,"^",fixed = T),
                                                                                      function(x){x[length(x)-1]})),"-"),function(x){x[2]})))
    ##################### SE
    ######## SE: upstream intron 300bp, whole exon, downsteam intron 300bp
    # 1: upstream intron
    bed6_upstream_intron <- data.frame(chr = diff_splicing_anno_SE_sel$chr,
                                       start = ifelse(diff_splicing_anno_SE_sel$strand == "+",
                                                      diff_splicing_anno_SE_sel$start-region,
                                                      diff_splicing_anno_SE_sel$start),
                                       end = ifelse(diff_splicing_anno_SE_sel$strand == "+",
                                                    diff_splicing_anno_SE_sel$start,
                                                    diff_splicing_anno_SE_sel$start+region),
                                       name = diff_splicing_anno_SE_sel$full_id,
                                       score = rep(1,nrow(diff_splicing_anno_SE_sel)),
                                       strand = diff_splicing_anno_SE_sel$strand)
    
    bed6_upstream_region <- data.frame(chr = diff_splicing_anno_SE_sel$chr,
                                       start = ifelse(diff_splicing_anno_SE_sel$strand == "+",
                                                      diff_splicing_anno_SE_sel$start-region,
                                                      diff_splicing_anno_SE_sel$start-region),
                                       end = ifelse(diff_splicing_anno_SE_sel$strand == "+",
                                                    diff_splicing_anno_SE_sel$start+region,
                                                    diff_splicing_anno_SE_sel$start+region),
                                       name = diff_splicing_anno_SE_sel$full_id,
                                       score = rep(1,nrow(diff_splicing_anno_SE_sel)),
                                       strand = diff_splicing_anno_SE_sel$strand)
    
    # 2: full exon region
    bed6_exon <- data.frame(chr = diff_splicing_anno_SE_sel$chr,
                            start = ifelse(diff_splicing_anno_SE_sel$strand == "+",
                                           diff_splicing_anno_SE_sel$start-region,
                                           diff_splicing_anno_SE_sel$end-region),
                            end = ifelse(diff_splicing_anno_SE_sel$strand == "+",
                                         diff_splicing_anno_SE_sel$end+region,
                                         diff_splicing_anno_SE_sel$start+region),
                            name = diff_splicing_anno_SE_sel$full_id,
                            score = rep(1,nrow(diff_splicing_anno_SE_sel)),
                            strand = diff_splicing_anno_SE_sel$strand)
    
    # 3: downstream intron
    bed6_downstream_intron <- data.frame(chr = diff_splicing_anno_SE_sel$chr,
                                         start = ifelse(diff_splicing_anno_SE_sel$strand == "+",
                                                        diff_splicing_anno_SE_sel$end,
                                                        diff_splicing_anno_SE_sel$end-region),
                                         end = ifelse(diff_splicing_anno_SE_sel$strand == "+",
                                                      diff_splicing_anno_SE_sel$end+region,
                                                      diff_splicing_anno_SE_sel$end),
                                         name = diff_splicing_anno_SE_sel$full_id,
                                         score = rep(1,nrow(diff_splicing_anno_SE_sel)),
                                         strand = diff_splicing_anno_SE_sel$strand)
    
    bed6_downstream_region <- data.frame(chr = diff_splicing_anno_SE_sel$chr,
                                         start = ifelse(diff_splicing_anno_SE_sel$strand == "+",
                                                        diff_splicing_anno_SE_sel$end-region,
                                                        diff_splicing_anno_SE_sel$end-region),
                                         end = ifelse(diff_splicing_anno_SE_sel$strand == "+",
                                                      diff_splicing_anno_SE_sel$end+region,
                                                      diff_splicing_anno_SE_sel$end+region),
                                         name = diff_splicing_anno_SE_sel$full_id,
                                         score = rep(1,nrow(diff_splicing_anno_SE_sel)),
                                         strand = diff_splicing_anno_SE_sel$strand)
    
    # full length in +/- region
    bed6 <- data.frame(chr = diff_splicing_anno_SE_sel$chr,
                       start = ifelse(diff_splicing_anno_SE_sel$strand == "+",
                                      diff_splicing_anno_SE_sel$start-region,
                                      diff_splicing_anno_SE_sel$end-region),
                       end = ifelse(diff_splicing_anno_SE_sel$strand == "+",
                                    diff_splicing_anno_SE_sel$end+region,
                                    diff_splicing_anno_SE_sel$start+region),
                       name = diff_splicing_anno_SE_sel$full_id,
                       score = rep(1,nrow(diff_splicing_anno_SE_sel)),
                       strand = diff_splicing_anno_SE_sel$strand)
    
    write.table(bed6_upstream_intron,paste(diff_age_deepbind_dir,tissue_dir,"/bed6_300bp_",choose_sex,"/SE_upstream_intron.bed6",sep = ""),sep = "\t",quote = F,col.names = F,row.names = F)
    write.table(bed6_exon,paste(diff_age_deepbind_dir,tissue_dir,"/bed6_300bp_",choose_sex,"/SE_exon.bed6",sep = ""),sep = "\t",quote = F,col.names = F,row.names = F)
    write.table(bed6_downstream_intron,paste(diff_age_deepbind_dir,tissue_dir,"/bed6_300bp_",choose_sex,"/SE_downstream_intron.bed6",sep = ""),sep = "\t",quote = F,col.names = F,row.names = F)
    write.table(bed6,paste(diff_age_deepbind_dir,tissue_dir,"/bed6_300bp_",choose_sex,"/SE_full_region.bed6",sep = ""),sep = "\t",quote = F,col.names = F,row.names = F)
    write.table(bed6_upstream_region,paste(diff_age_deepbind_dir,tissue_dir,"/bed6_300bp_",choose_sex,"/SE_upstream_region.bed6",sep = ""),sep = "\t",quote = F,col.names = F,row.names = F)
    write.table(bed6_downstream_region,paste(diff_age_deepbind_dir,tissue_dir,"/bed6_300bp_",choose_sex,"/SE_downstream_region.bed6",sep = ""),sep = "\t",quote = F,col.names = F,row.names = F)
    
  }
  
  if(nrow(diff_splicing_anno_RI_sel)>0){
    diff_splicing_anno_RI_sel$splice_out <- unlist(lapply(strsplit(diff_splicing_anno_RI_sel$full_id,"@"),function(x){x[6]}))
    diff_splicing_anno_RI_sel$start <-  as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(diff_splicing_anno_RI_sel$splice_out,"^",fixed = T),
                                                                                        function(x){x[2]})),"-"),function(x){x[1]})))
    diff_splicing_anno_RI_sel$end <-  as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(diff_splicing_anno_RI_sel$splice_out,"^",fixed = T),
                                                                                      function(x){x[length(x)-1]})),"-"),function(x){x[2]})))
    
    ##################### RI
    ####### RI: upstream exon 300bp, upstream intron 200 bp; downstream intron 300bp, downstream exon 300bp
    # 1: upstream exon
    bed6_upstream_exon <- data.frame(chr = diff_splicing_anno_RI_sel$chr,
                                     start = ifelse(diff_splicing_anno_RI_sel$strand == "+",
                                                    diff_splicing_anno_RI_sel$start-region,
                                                    diff_splicing_anno_RI_sel$start),
                                     end = ifelse(diff_splicing_anno_RI_sel$strand == "+",
                                                  diff_splicing_anno_RI_sel$start,
                                                  diff_splicing_anno_RI_sel$start+region),
                                     name = diff_splicing_anno_RI_sel$full_id,
                                     score = rep(1,nrow(diff_splicing_anno_RI_sel)),
                                     strand = diff_splicing_anno_RI_sel$strand)
    # 2: upstream intron
    bed6_upstream_intron <- data.frame(chr = diff_splicing_anno_RI_sel$chr,
                                       start = ifelse(diff_splicing_anno_RI_sel$strand == "+",
                                                      diff_splicing_anno_RI_sel$start,
                                                      diff_splicing_anno_RI_sel$start-region),
                                       end = ifelse(diff_splicing_anno_RI_sel$strand == "+",
                                                    diff_splicing_anno_RI_sel$start+region,
                                                    diff_splicing_anno_RI_sel$start),
                                       name = diff_splicing_anno_RI_sel$full_id,
                                       score = rep(1,nrow(diff_splicing_anno_RI_sel)),
                                       strand = diff_splicing_anno_RI_sel$strand)
    # 3: downstream intron
    bed6_downstream_intron <- data.frame(chr = diff_splicing_anno_RI_sel$chr,
                                         start = ifelse(diff_splicing_anno_RI_sel$strand == "+",
                                                        diff_splicing_anno_RI_sel$end-region,
                                                        diff_splicing_anno_RI_sel$end),
                                         end = ifelse(diff_splicing_anno_RI_sel$strand == "+",
                                                      diff_splicing_anno_RI_sel$end,
                                                      diff_splicing_anno_RI_sel$end+region),
                                         name = diff_splicing_anno_RI_sel$full_id,
                                         score = rep(1,nrow(diff_splicing_anno_RI_sel)),
                                         strand = diff_splicing_anno_RI_sel$strand)
    # 4: downstream exon
    bed6_downstream_exon <- data.frame(chr = diff_splicing_anno_RI_sel$chr,
                                       start = ifelse(diff_splicing_anno_RI_sel$strand == "+",
                                                      diff_splicing_anno_RI_sel$end,
                                                      diff_splicing_anno_RI_sel$end-region),
                                       end = ifelse(diff_splicing_anno_RI_sel$strand == "+",
                                                    diff_splicing_anno_RI_sel$end+region,
                                                    diff_splicing_anno_RI_sel$end),
                                       name = diff_splicing_anno_RI_sel$full_id,
                                       score = rep(1,nrow(diff_splicing_anno_RI_sel)),
                                       strand = diff_splicing_anno_RI_sel$strand)
    # upstream region
    bed6_upstream <- data.frame(chr = diff_splicing_anno_RI_sel$chr,
                                start = ifelse(diff_splicing_anno_RI_sel$strand == "+",
                                               diff_splicing_anno_RI_sel$start-region,
                                               diff_splicing_anno_RI_sel$start-region),
                                end = ifelse(diff_splicing_anno_RI_sel$strand == "+",
                                             diff_splicing_anno_RI_sel$start+region,
                                             diff_splicing_anno_RI_sel$start+region),
                                name = diff_splicing_anno_RI_sel$full_id,
                                score = rep(1,nrow(diff_splicing_anno_RI_sel)),
                                strand = diff_splicing_anno_RI_sel$strand)
    # downstream region
    bed6_downstream <- data.frame(chr = diff_splicing_anno_RI_sel$chr,
                                  start = ifelse(diff_splicing_anno_RI_sel$strand == "+",
                                                 diff_splicing_anno_RI_sel$end-region,
                                                 diff_splicing_anno_RI_sel$end-region),
                                  end = ifelse(diff_splicing_anno_RI_sel$strand == "+",
                                               diff_splicing_anno_RI_sel$end+region,
                                               diff_splicing_anno_RI_sel$end+region),
                                  name = diff_splicing_anno_RI_sel$full_id,
                                  score = rep(1,nrow(diff_splicing_anno_RI_sel)),
                                  strand = diff_splicing_anno_RI_sel$strand)
    
    write.table(bed6_upstream_exon,paste(diff_age_deepbind_dir,tissue_dir,"/bed6_300bp_",choose_sex,"/RI_upstream_exon.bed6",sep = ""),sep = "\t",quote = F,col.names = F,row.names = F)
    write.table(bed6_upstream_intron,paste(diff_age_deepbind_dir,tissue_dir,"/bed6_300bp_",choose_sex,"/RI_upstream_intron.bed6",sep = ""),sep = "\t",quote = F,col.names = F,row.names = F)
    write.table(bed6_downstream_intron,paste(diff_age_deepbind_dir,tissue_dir,"/bed6_300bp_",choose_sex,"/RI_downstream_intron.bed6",sep = ""),sep = "\t",quote = F,col.names = F,row.names = F)
    write.table(bed6_downstream_exon,paste(diff_age_deepbind_dir,tissue_dir,"/bed6_300bp_",choose_sex,"/RI_downstream_exon.bed6",sep = ""),sep = "\t",quote = F,col.names = F,row.names = F)
    write.table(bed6_upstream,paste(diff_age_deepbind_dir,tissue_dir,"/bed6_300bp_",choose_sex,"/RI_upstream_region.bed6",sep = ""),sep = "\t",quote = F,col.names = F,row.names = F)
    write.table(bed6_downstream,paste(diff_age_deepbind_dir,tissue_dir,"/bed6_300bp_",choose_sex,"/RI_downstream_region.bed6",sep = ""),sep = "\t",quote = F,col.names = F,row.names = F)
    
  }
  
  if(nrow(diff_splicing_anno_A5SS_sel)>0){
    diff_splicing_anno_A5SS_sel$start <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(diff_splicing_anno_A5SS_sel$full_id,"@",fixed = T),
                                                                                         function(x){x[5]})),"^",fixed = T),function(x){x[1]})))
    diff_splicing_anno_A5SS_sel$end <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(diff_splicing_anno_A5SS_sel$full_id,"@",fixed = T),
                                                                                       function(x){x[6]})),"^",fixed = T),function(x){x[1]})))
    
    ##################### A5SS
    # A5SS: 1st ss upsteam exon 300bp, 1st ss to 2nd ss, 2nd ss downstream intron 300bp
    # full region
    bed6 <- data.frame(chr = diff_splicing_anno_A5SS_sel$chr,
                       start = ifelse(diff_splicing_anno_A5SS_sel$strand == "+",
                                      diff_splicing_anno_A5SS_sel$start-region,
                                      diff_splicing_anno_A5SS_sel$end-region),
                       end = ifelse(diff_splicing_anno_A5SS_sel$strand == "+",
                                    diff_splicing_anno_A5SS_sel$end+region,
                                    diff_splicing_anno_A5SS_sel$start+region),
                       name = diff_splicing_anno_A5SS_sel$full_id,
                       score = rep(1,nrow(diff_splicing_anno_A5SS_sel)),
                       strand = diff_splicing_anno_A5SS_sel$strand)
    write.table(bed6,paste(diff_age_deepbind_dir,tissue_dir,"/bed6_300bp_",choose_sex,"/A5SS_full_region.bed6",sep = ""),sep = "\t",quote = F,col.names = F,row.names = F)
    
  }
  
  if(nrow(diff_splicing_anno_A3SS_sel)>0){
    diff_splicing_anno_A3SS_sel$start <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(diff_splicing_anno_A3SS_sel$full_id,"@",fixed = T),
                                                                                                                function(x){x[6]})),"^",fixed = T),function(x){x[2]})),"-"),
                                                                  function(x){x[1]})))
    diff_splicing_anno_A3SS_sel$end <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(diff_splicing_anno_A3SS_sel$full_id,"@",fixed = T),
                                                                                                              function(x){x[5]})),"^",fixed = T),function(x){x[2]})),"-"),
                                                                function(x){x[1]})))
    ##################### A3SS
    # A3SS: 1st ss upsteam intron 300bp, 1st ss to 2nd ss, 2nd ss downstream exon 300bp
    # full region
    bed6 <- data.frame(chr = diff_splicing_anno_A3SS_sel$chr,
                       start = ifelse(diff_splicing_anno_A3SS_sel$strand == "+",
                                      diff_splicing_anno_A3SS_sel$start-region,
                                      diff_splicing_anno_A3SS_sel$end-region),
                       end = ifelse(diff_splicing_anno_A3SS_sel$strand == "+",
                                    diff_splicing_anno_A3SS_sel$end+region,
                                    diff_splicing_anno_A3SS_sel$start+region),
                       name = diff_splicing_anno_A3SS_sel$full_id,
                       score = rep(1,nrow(diff_splicing_anno_A3SS_sel)),
                       strand = diff_splicing_anno_A3SS_sel$strand)
    write.table(bed6,paste(diff_age_deepbind_dir,tissue_dir,"/bed6_300bp_",choose_sex,"/A3SS_full_region.bed6",sep = ""),sep = "\t",quote = F,col.names = F,row.names = F)
  }
}



######################################################################
##### run bedtools
# conda activate base
### shell run bedtools
# bash run_bedtools_batch.sh
#bedtools getfasta -fi GRCh38.p13.genome.fa -name -s -bed $bed6 -fo $fa ;
######################################################################



######################################################################
# make files to deepbind
######################################################################
# deepbind subsequences window <- 40
# direction: 5' -> 3'
library("stringr")

split_seq_window <- function(x){
  len <- str_length(x)
  seqs <- c()
  for (i in 1:(len-window+1)) {
    seq <- substr(x,i,i+window-1)
    seqs <-c(seqs,seq)
  }
  return(seqs)
}
tissue_dir <- ""
print(tissue_dir)
for (choose_sex in c("female","male")) {
  print(choose_sex)
  dirs <- gsub(".fa","",list.files(paste(diff_age_deepbind_dir,tissue_dir,"/fa_300bp_",choose_sex,"/",sep = "")))
  for (dir in dirs) {
    print(dir)
    
    system(paste("mkdir ",diff_age_deepbind_dir,gsub("\\)","\\\\)",gsub("\\(","\\\\(",tissue_dir)),"/seq_window_40bp_",choose_sex,"/",dir,sep = ""))
    if((file.info(paste(diff_age_deepbind_dir,tissue_dir,"/fa_300bp_",choose_sex,"/",dir,".fa",sep = ""))$size <= 2)){next}
    
    fa <- read.delim(paste(diff_age_deepbind_dir,tissue_dir,"/fa_300bp_",choose_sex,"/",dir,".fa",sep = ""),header = F)
    bed6 <- read.delim(paste(diff_age_deepbind_dir,tissue_dir,"/bed6_300bp_",choose_sex,"/",dir,".bed6",sep = ""),header = F)
    df_fa <- data.frame(seqs = fa[seq(0,nrow(fa),2),],
                        strand = bed6$V6,
                        events = bed6$V4)
    events_seqs <- mclapply(df_fa$seqs,split_seq_window,mc.cores = 2)
    names(events_seqs) <- df_fa$events
    
    for (i in names(events_seqs)) {
      write.table(as.matrix(events_seqs[[i]]),paste(diff_age_deepbind_dir,tissue_dir,"/seq_window_40bp_",choose_sex,"/",dir,"/",dir,"-",i,".seq",sep = ""),
                  quote = F,col.names = F,row.names = F)
    }
  }
}

######################################################################
#### run deepbind 

# bash run_deepbind.sh
# for j in `ls ./${tissue}/seq_window_40bp_male/${i}/${i}*`;
# do
# events=${j/.seq/};
# events=${events/\.\/${tissue}\/seq_window_40bp_male\/${i}\/${i}\-/};
# ids=./${tissue}/${tissue}_deepbind-rbp-hs_male.ids;
# db=./${tissue}/deepbind_diff_exp_rbp_40bp_male/${i}/${i}-${events}-deepbind.txt;
# deepbind $ids $j > $db ;
######################################################################

