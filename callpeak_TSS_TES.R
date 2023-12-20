options(stringsAsFactors = F)
options(scipen = 200)
library(data.table)
library(stringr)
library(dplyr)
library('BSgenome.Hsapiens.UCSC.hg38')
library(bedtoolsr)
library(CAGEr)
args=commandArgs(T)


f<-fread(args[1],sep = '\t',fill=T) %>% as.data.frame()
f$umi_map_reads <- paste0(f$UMI,'_',f$V6,'__',f$V13)
f$map_reads <- paste0(f$V6,'__',f$V13)
f_dedup1 <- f[!duplicated(f$umi_map_reads),]
f_dedup2 <- f[!duplicated(f$map_reads),]
f <- f_dedup1[,-c(22,23)]

UMI_TSS <- unique(f$UMI[grep('TSS',f$V1)])
UMI_TES <- unique(f$UMI[grep('TES',f$V1)])
UMI_exon <- unique(f$UMI[grep('exon',f$V1)])
UMI_exon_TSS <- unique(f$UMI[grep('exonHead',f$V1)])
UMI_exon_TES <- unique(f$UMI[grep('exonTail',f$V1)])
UMI_complete<-intersect(UMI_TSS,intersect(UMI_TES,UMI_exon)) %>% unique()


f$strand <- '+'
f$strand[grep('16',f$V2)]<-'-'
f_plus<-f[grep('0',f$V2),]
f_plus_TSS<-f_plus[grep('TSS',f_plus$V1),]
f_plus_TES<-f_plus[grep('TES',f_plus$V1),]
f_plus_TSS$TSS<-str_remove(f_plus_TSS$V13,'-.*$')
f_plus_TES$TES<-str_remove(f_plus_TES$V13,'^.*-')
f_minus<-f[grep('16',f$V2),]
f_minus_TSS<-f_minus[grep('TSS',f_minus$V1),]
f_minus_TES<-f_minus[grep('TES',f_minus$V1),]
f_minus_TSS$TSS<-str_remove(f_minus_TSS$V13,'^.*-')
f_minus_TES$TES<-str_remove(f_minus_TES$V13,'-.*$')
f_tss<-rbind(f_plus_TSS,f_minus_TSS)
f_tes<-rbind(f_plus_TES,f_minus_TES)
f_tss$TSS_strand<-paste0(f_tss$TSS,'_',f_tss$strand)
f_tes$TES_strand<-paste0(f_tes$TES,'_',f_tes$strand)

library_UMI <- data.frame(UMI=unique(f$UMI),TSS=0,exon=0,TES=0)
library_UMI$TSS[library_UMI$UMI%in% UMI_TSS] <- 1
library_UMI$TES[library_UMI$UMI%in% UMI_TES] <- 1
library_UMI$exon[library_UMI$UMI%in% UMI_exon] <- 1


if(!dir.exists(paste0(dir,'temp'))){
  dir.create(paste0(dir,'temp'))
}
tss_bed6<-data.frame(f_tss$V3,round(as.numeric(f_tss$TSS)-1),as.numeric(f_tss$TSS),f_tss$UMI,'0',f_tss$strand)
tes_bed6<-data.frame(f_tes$V3,round(as.numeric(f_tes$TES)-1),as.numeric(f_tes$TES),f_tes$UMI,'0',f_tes$strand)
write.table(tss_bed6,file = paste0(dir,'temp','/tss_bed6.bed'),quote = F,row.names = F,col.names = F,sep = "\t")
write.table(tes_bed6,file = paste0(dir,'temp','/tes_bed6.bed'),quote = F,row.names = F,col.names = F,sep = "\t")

myCAGEset <- new("CAGEset", genomeName = genomeName,
                 inputFiles = paste0(dir,'temp','/tss_bed6.bed'),
                 inputFilesType = "bed",sampleLabels = "sample")
getCTSS(myCAGEset,removeFirstG = FALSE,correctSystematicG = FALSE, nrCores = 10)
normalizeTagCount(myCAGEset,method = "none",T = 10^6)
clusterCTSS(object = myCAGEset, method = "distclu", threshold = 0, nrPassThreshold = 1, thresholdIsTpm = TRUE,
            removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = T, nrCores = 40)
tc_tss <- tagClusters(myCAGEset)[[1]]

df_cluster_tss<-NULL
for (i in 1:nrow(tc_tss)) {
  df<-tc_tss[i,]
  end_loc<-(df$start+1):df$end
  df_cluster<-cbind(end_loc,paste0('cluster',df$cluster),as.character(df$strand)) %>% as.data.frame()
  colnames(df_cluster)<-c('end_loc','tss_cluster','strand')
  df_cluster$range <- paste0((df$start+1),'-',df$end)
  df_cluster$headmost_tss <- ifelse(df$strand=='+',df$start+1,df$end)
  
  df_cluster_tss<-rbind(df_cluster_tss,df_cluster)
}
df_cluster_tss$TSS_strand<-paste0(df_cluster_tss$end_loc,'_',df_cluster_tss$strand)
cluster_tss<-merge(f_tss,df_cluster_tss[,c(2,4,5,6)],by='TSS_strand',all.x=T,all.y=F)
cluster_tss<-cluster_tss[,-1]
cluster_tss$UMI_cluster<-paste0(cluster_tss$UMI,'-',cluster_tss$tss_cluster)

myCAGEset <- new("CAGEset", genomeName = genomeName,
                 inputFiles = paste0(dir,'temp','/tes_bed6.bed'),
                 inputFilesType = "bed",sampleLabels = "sample")
getCTSS(myCAGEset,removeFirstG = FALSE,correctSystematicG = FALSE, nrCores = 10)
normalizeTagCount(myCAGEset,method = "none",T = 10^6)
clusterCTSS(object = myCAGEset, method = "distclu", threshold = 0, nrPassThreshold = 1, thresholdIsTpm = TRUE,
            removeSingletons = FALSE,keepSingletonsAbove = 1, maxDist = 20,useMulticore = T, nrCores = 40)
tc_tes <- tagClusters(myCAGEset)[[1]]

df_cluster_tes<-NULL
for (i in 1:nrow(tc_tes)) {
  df<-tc_tes[i,]
  end_loc<-(df$start+1):df$end
  df_cluster<-cbind(end_loc,paste0('cluster',df$cluster),as.character(df$strand)) %>% as.data.frame()
  colnames(df_cluster)<-c('end_loc','tes_cluster','strand')
  df_cluster$range <- paste0((df$start+1),'-',df$end)
  df_cluster$hindmost_tes <- ifelse(df$strand=='+',df$end,df$start+1)
  
  df_cluster_tes<-rbind(df_cluster_tes,df_cluster)
}
df_cluster_tes$TES_strand<-paste0(df_cluster_tes$end_loc,'_',df_cluster_tes$strand)
cluster_tes<-merge(f_tes,df_cluster_tes[,c(2,4,5,6)],by='TES_strand',all.x=T,all.y=F)
cluster_tes<-cluster_tes[,-1]
cluster_tes$UMI_cluster<-paste0(cluster_tes$UMI,'-',cluster_tes$tes_cluster)


table_tss<-table(cluster_tss$UMI_cluster) %>% sort(decreasing = T)
df_tss_fre <- as.data.frame(table_tss)
colnames(df_tss_fre) <- c('UMI_cluster','count')
df_tss_fre$UMI <- str_remove(df_tss_fre$UMI_cluster,'-cluster.*$')
df_tss_fre_dedup <- df_tss_fre[!duplicated(df_tss_fre$UMI),]
df_tss_fre_dedup <- merge(df_tss_fre_dedup,unique(cluster_tss[,25:27]),by='UMI_cluster',all.x=T,all.y=F)
multi_tss_umi<-df_tss_fre$UMI[duplicated(df_tss_fre$UMI)] %>% unique()


table_tes<-table(cluster_tes$UMI_cluster) %>% sort(decreasing = T)
df_tes_fre <- as.data.frame(table_tes)
colnames(df_tes_fre) <- c('UMI_cluster','count')
df_tes_fre$UMI <- str_remove(df_tes_fre$UMI_cluster,'-cluster.*$')
df_tes_fre_dedup <- df_tes_fre[!duplicated(df_tes_fre$UMI),]
df_tes_fre_dedup <- merge(df_tes_fre_dedup,unique(cluster_tes[,25:27]),by='UMI_cluster',all.x=T,all.y=F)
multi_tes_umi<-df_tes_fre$UMI[duplicated(df_tes_fre$UMI)] %>% unique()
multi_umi<-union(multi_tss_umi,multi_tes_umi)
conflict_UMI <- multi_umi

library_UMI$stage <- 'correct'
library_UMI$stage[library_UMI$UMI %in% conflict_UMI] <- 'conflict'
library_UMI <- merge(library_UMI,df_tss_fre_dedup[,c(3,5)],by='UMI',all.x=T,all.y=F)
library_UMI <- merge(library_UMI,df_tes_fre_dedup[,c(3,5)],by='UMI',all.x=T,all.y=F)
library_UMI$reads_count <- count[library_UMI$UMI]
save(library_UMI,file='library_UMI.Rdata')




