options(stringsAsFactors = F)
library(dplyr)
library(stringr)
library(bedtoolsr)
library(data.table)

file <- args[1]
chr <- strsplit(file,'/') %>% unlist()
chr <- chr[length(chr)-1]
dir_ref <- str_remove(file,'assigned_isoforms.*$')
if(!dir.exists(paste0(dir_ref,'assigned_isoforms/',chr,'/bed'))){
  dir.create(paste0(dir_ref,'assigned_isoforms/',chr,'/bed'))
}

gtf <-  read.table("index/hg38_SIRVercc/hg38_SIRVercc.gtf",sep = "\t")
gtf_tran <- gtf[gtf$V3=='transcript',]
genes <- str_remove(gtf_tran$V9,'gene_id ') %>% str_remove(';.*$')  
transcripts <- str_remove(gtf_tran$V9,'^.*transcript_id ') %>% str_remove(';.*$') 
names(genes) <- transcripts
gtf_gene <- gtf[gtf$V3=='gene',]
gtf_gene$V9 <- str_remove(gtf_gene$V9,'gene_id ') %>% str_remove(';.*$')

f_all <- fread(file,sep='\t',fill=T,header = F) %>% as.data.frame()
f_all$all_UMI <- paste0(f_all$V2,'____',f_all$V5)
f_all$UMI <- paste0(f_all$V2,'____',genes[str_remove(f_all$V5,',.*$')])

load(paste0(str_remove(file,'assign.*$'),'keptReads/',chr,'/library_UMI.Rdata'))
library_UMI <- library_UMI[library_UMI$reads_count>=1,]
library_UMI <- library_UMI[,-8]
f_all <- merge(f_all,library_UMI[library_UMI$stage=='correct',],by.y='UMI',all.x=F,all.y=F)
f_all$all_UMI <- paste0(f_all$all_UMI,'_TSS',f_all$TSS,'_exon',f_all$exon,'_TES',f_all$TES,'_TSSloc',f_all$TSSloc,'_TESloc',f_all$TESloc)
f_all <- f_all[,c(2,3,4,5,6,7,8,9,1)]
f_all$all_UMI <- paste0(f_all$V1,'--',f_all$all_UMI)

start<-Sys.time()
print(paste0('start time ',start))
cat('\n')

f_new <- data.frame()        
f_exon <- data.frame()       
f_gap <- data.frame()        
for (i in 1:nrow(f_all)) {
  
  cat(i,' in ',nrow(f_all),' UMIs\n')
  
  df <- f_all[i,]
  reads <- strsplit(df[,4],';') %>% unlist() 
  num <- strsplit(df[,4],';') %>% unlist() %>% strsplit(',') %>% unlist() %>% str_split('-') %>% unlist()  %>% as.numeric()
  min <- num %>% min()
  max <- num  %>% max()

  df_gene <- gtf_gene[gtf_gene$V9==str_remove(df$UMI,'^.*____'),]
  min <- ifelse(min+500<df_gene[,4],df_gene[,4],min)
  max <- ifelse(max-500>df_gene[,5],df_gene[,5],max)
  
  df_all <- data.frame(V1=df[,8],V2=num[seq(1,length(num),2)]-1,V3=num[seq(2,length(num),2)],V4=df[,8],V5='.',V6=df_gene[,7]) %>% unique()
  df_all <- df_all[order(df_all$V2),]  
  reads1 <- grep(',',reads,value = T)
  if(length(reads1)>0){
    gap <- lapply(reads1,function(x){grep(',',unlist(strsplit(x,'-')),value = T) %>% cbind()})
    gap1 <- do.call(rbind,gap) %>% str_split(',',simplify = T) %>% as.data.frame()
    gap2 <- gap1[order(gap1[,1],decreasing = F),] %>% unique()
    gap2 <- gap2[gap2$V1>min|gap2$V2<max,]
    gap2$V1 <- as.numeric(gap2$V1) ; gap2$V2 <- as.numeric(gap2$V2)
    gap2 <- gap2[gap2$V2-gap2$V1>10,]  
    if(nrow(gap2)>0){
      df_gap <- data.frame(V1=df[,8],V2=gap2[,1],V3=gap2[,2],V4=df[,8],V5='.',V6=df_gene[,7])             
      df_new <- data.frame(V1=df[,8],V2=c(min,gap2$V2)-1,V3=c(gap2$V1,max),V4=df[,8],V5='.',V6=df_gene[,7])
      f_new <- rbind(f_new,df_new)
      f_gap <- rbind(f_gap,df_gap)
      
    }else{
      df_new <- data.frame(V1=df[,8],V2=min-1,V3=max,V4=df[,8],V5='.',V6=df_gene[,7])  
      f_new <- rbind(f_new,df_new)
    }
    
    
  }else{
    df_new <- data.frame(V1=df[,8],V2=min-1,V3=max,V4=df[,8],V5='.',V6=df_gene[,7]) 
    f_new <- rbind(f_new,df_new)
  }
  
  f_exon <- rbind(f_exon,df_all)
  
  if(i%%5000==0) { 
    gc() }
  
  
}


f_gap1 <- f_gap
f_gap1$V2 <- f_gap1$V2+3
f_gap1$V3 <- f_gap1$V3-3
conflict <- bt.intersect(a=f_exon,b=f_gap1,s=T,wa=T,wb=T)

f_new <- f_new[! f_new$V4 %in% conflict$V4,]
end<-Sys.time()
print(paste0('end time ',end))
cat('\n')
print(paste0('used time ',end-start))
cat('\n')

f_new[,1] <- str_remove(f_new[,1],'_TSSloc.*$')
f_new[,4] <- str_remove(f_new[,4],'_TSSloc.*$')
f_new <- f_new[order(f_new[,3],decreasing = F),]
f_new <- f_new[order(f_new[,2],decreasing = F),] 
f_new <- f_new[order(f_new[,1],decreasing = F),] 
f_new <- f_new[f_new[,3]>f_new[,2],]

write.table(f_new,file = paste0(dir_ref,'assigned_isoforms/',chr,"/bed/for_merge.bed"),quote = F,row.names = F,col.names = F,sep = "\t")
system(paste0("cgat bed2bed --method=merge --merge-by-name -I ",dir_ref,'assigned_isoforms/',chr,"/bed/for_merge.bed"," -S ",dir_ref,'assigned_isoforms/',chr,"/bed/no_redundance.bed"))
system(paste0('python3 bed6Tobed12.py ',dir_ref,'assigned_isoforms/',chr,'/bed/no_redundance.bed > ',dir_ref,'assigned_isoforms/',chr,'/bed/all_exon_12.bed'))












