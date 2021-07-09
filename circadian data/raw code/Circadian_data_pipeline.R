rm(list = ls())
library(ggplot2)
library(pheatmap)
require(reshape2)
library(MetaCycle)
source("../../materials//scTFA_calculation.R")

themo_demo=theme(
  text = element_text(size=22),
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank(),
  axis.text=element_text(color='black'),
  plot.title = element_text(hjust = 0.5),
  title=element_text(size = 30),
  axis.title.x = element_text(size=24),
  axis.title.y = element_text(size=24),
  axis.text.x = element_text(size=22,angle = 45,hjust = 1),
  axis.text.y = element_text(size=22))

data_dir="../raw data/"
result_dir="../results/"
date_analysis="analyzed20210201_2(metacycle)"
dir.create(paste(result_dir))

#import RNA data
RNA_data=read.table(paste(data_dir,"GSE73554_WT_AL_Intron_Exon_RFP.txt",sep=""),header=T)
RNA_RPKM=RNA_data
RNA_RPKM[,3:ncol(RNA_RPKM)]=2^RNA_RPKM[,3:ncol(RNA_RPKM)]#original data is log2(RPKM)

#aggregate rows with the same gene symbol by summing
RNA_RPKM_aggregated=aggregate(RNA_RPKM[,3:ncol(RNA_RPKM)],by=list(Gene_Symbol=RNA_RPKM$Gene_Symbol),FUN = sum)
rownames(RNA_RPKM_aggregated)=RNA_RPKM_aggregated$Gene_Symbol
RNA_RPKM_aggregated=RNA_RPKM_aggregated[,-1]
# colnames(RNA_RPKM_aggregated)

#split by time points and intron/exon
condition="WT_AL"
type=c("Intron","Exon")
time=c("00","02","04","06","08","10","12","14","16","18","20","22")

RNA_RPKM_split=matrix(,nrow(RNA_RPKM_aggregated),length(type)*length(time),
                      dimnames = list(rownames(RNA_RPKM_aggregated),paste(rep(type,each=length(time)),time,sep="_")))
for(type_i in type){
  for(time_i in time){
    RNA_RPKM_split[,paste(type_i,time_i,sep="_")]=apply(RNA_RPKM_aggregated[,grep(paste(condition,type_i,time_i,sep="_"),
                                                                                  colnames(RNA_RPKM_aggregated),value=T)],
                                                        1,mean)
  }
}
write.csv(RNA_RPKM_split,paste(result_dir,"RNA_RPKM_split",date_analysis,".csv",sep = ""))

######## alternative calculation
if(F){
  # use only one type (A,B,C,D) to calculate
  subtypes = c("A","B","C","D")
  subtype = subtypes[1]
  RNA_RPKM_split=matrix(,nrow(RNA_RPKM_aggregated),length(type)*length(time),
                        dimnames = list(rownames(RNA_RPKM_aggregated),paste(rep(type,each=length(time)),time,sep="_")))
  for(type_i in type){
    for(time_i in time){
      RNA_RPKM_split[,paste(type_i,time_i,sep="_")]=RNA_RPKM_aggregated[,grep(paste(condition,type_i,time_i,subtype,sep="_"),
                                                                              colnames(RNA_RPKM_aggregated),value=T)]
    }
  }
  write.csv(RNA_RPKM_split,paste(result_dir,"RNA_RPKM_split",date_analysis,".csv",sep = ""))
  
}
if(F){
  # remove one subtype
  cir_TF_list = NULL
  for(subtype_i in 1:4){
    
    subtypes = c("A","B","C","D")
    subtype = subtypes[subtype_i]
    RNA_RPKM_split=matrix(,nrow(RNA_RPKM_aggregated),length(type)*length(time),
                          dimnames = list(rownames(RNA_RPKM_aggregated),paste(rep(type,each=length(time)),time,sep="_")))
    for(type_i in type){
      for(time_i in time){
        target_types = grep(paste(condition,type_i,time_i,sep="_"),colnames(RNA_RPKM_aggregated),value=T)
        target_types = setdiff(target_types,grep(paste(condition,type_i,time_i,subtype,sep="_"),colnames(RNA_RPKM_aggregated),value=T))
        RNA_RPKM_split[,paste(type_i,time_i,sep="_")]=apply(RNA_RPKM_aggregated[,target_types],1,mean)
        
      }
    }
    write.csv(RNA_RPKM_split,paste(result_dir,"RNA_RPKM_split",date_analysis,".csv",sep = ""))
    
    
    if(T){
      
      # import GRN data
      dorothea_logic_matrix=read.csv("../../materials/mouse_GRN_classAB_activation_only/dorothea_mouse_AB.csv",
                                     row.names = 1, check.names = F)
      dorothea_logic_matrix=as.matrix(dorothea_logic_matrix)
      
      # TFA of RNA
      GA_RNA=RNA_RPKM_split
      GA_RNA=GA_RNA[rowSums(GA_RNA)>0,]
      
      TFA_RNA=scTFA_calculation(as.matrix(GA_RNA),dorothea_logic_matrix,zscore=F,
                                gene_weight_method=NULL)
      write.csv(TFA_RNA,paste(result_dir,"TFA_RNA data",date_analysis,".csv",sep = ""))
      
      TFA_RNA_Intron=TFA_RNA[,grep("Intron",colnames(TFA_RNA),value=T)]
      write.csv(TFA_RNA_Intron,paste(result_dir,"TFA_RNA_Intron",date_analysis,".csv",sep = ""))
      
      TFA_RNA_Exon=TFA_RNA[,grep("Exon",colnames(TFA_RNA),value=T)]
      write.csv(TFA_RNA_Exon,paste(result_dir,"TFA_RNA_Exon",date_analysis,".csv",sep = ""))
      
      #calculate circadian related parameters with MetaCycle
      Intron_meta2d_result=meta2d(infile=paste(result_dir,"TFA_RNA_Intron",date_analysis,".csv",sep = ""), filestyle="csv", 
                                  outdir=result_dir,
                                  timepoints=as.numeric(sub("Intron_","",colnames(TFA_RNA_Intron))),
                                  minper = 24,maxper = 24,ARSdefaultPer=24,
                                  outputFile=F)
      Intron_meta2d_result=Intron_meta2d_result$meta
      
      Exon_meta2d_result=meta2d(infile=paste(result_dir,"TFA_RNA_Exon",date_analysis,".csv",sep = ""), filestyle="csv", result_dir,
                                timepoints=as.numeric(sub("Exon_","",colnames(TFA_RNA_Exon))),
                                minper = 24,maxper = 24,ARSdefaultPer=24,
                                outputFile=F)
      Exon_meta2d_result=Exon_meta2d_result$meta
      
      Intron_Exon_meta2d_result=merge(Intron_meta2d_result,Exon_meta2d_result,by.x="CycID",by.y="CycID",all=F,
                                      suffixes=c("_Intron","_Exon"))
      write.csv(Intron_Exon_meta2d_result,paste(result_dir,"Intron_Exon_meta2d_result",date_analysis,".csv",sep =""))
      
      #choose TFs with circadian TFA from intron and exon (ARS p<0.05)
      Intron_Exon_meta2d_result_subset=Intron_Exon_meta2d_result[(Intron_Exon_meta2d_result$meta2d_pvalue_Intron<0.05)&
                                                                   (Intron_Exon_meta2d_result$meta2d_pvalue_Exon<0.05),]
      Intron_Exon_circadian_TF_vector=Intron_Exon_meta2d_result_subset$CycID
      
      # record the target number of the TFs
      # overlap between dorothea targets and GA_RNA data
      overlap_targets = intersect(rownames(GA_RNA),colnames(dorothea_logic_matrix))
      n_targets = rowSums(dorothea_logic_matrix[,overlap_targets])[Intron_Exon_circadian_TF_vector]
      
      phase_plot_data_frame=Intron_Exon_meta2d_result_subset[,c("CycID","meta2d_phase_Intron","meta2d_phase_Exon"),]
      phase_plot_data_frame = cbind(phase_plot_data_frame,n_targets)
      phase_plot_data_frame$Exon_minus_Intron_phase=phase_plot_data_frame$meta2d_phase_Exon-phase_plot_data_frame$meta2d_phase_Intron 
      
      phase_plot_data_frame$Exon_minus_Intron_phase[phase_plot_data_frame$Exon_minus_Intron_phase<(-12)]=
        phase_plot_data_frame$Exon_minus_Intron_phase[phase_plot_data_frame$Exon_minus_Intron_phase<(-12)]+24
      phase_plot_data_frame$Exon_minus_Intron_phase[phase_plot_data_frame$Exon_minus_Intron_phase>12]=
        phase_plot_data_frame$Exon_minus_Intron_phase[phase_plot_data_frame$Exon_minus_Intron_phase>12]-24
      
      phase_plot_data_frame=phase_plot_data_frame[order(phase_plot_data_frame$Exon_minus_Intron_phase,decreasing = T),]
      phase_plot_data_frame$CycID=factor(phase_plot_data_frame$CycID,levels=phase_plot_data_frame$CycID)
      
    }
    
    cir_TF_list[[subtype_i]] = Intron_Exon_circadian_TF_vector
    print(subtype_i)
  }
  robustness = as.data.frame(sort(table(unlist(cir_TF_list)),decreasing = T))
  
}

if(T){
  
  # import GRN data
  dorothea_logic_matrix=read.csv("../../materials/mouse_GRN_classAB_activation_only/dorothea_mouse_AB.csv",
                                 row.names = 1, check.names = F)
  dorothea_logic_matrix=as.matrix(dorothea_logic_matrix)
  
  # TFA of RNA
  GA_RNA=RNA_RPKM_split
  GA_RNA=GA_RNA[rowSums(GA_RNA)>0,]
  
  TFA_RNA=scTFA_calculation(as.matrix(GA_RNA),dorothea_logic_matrix,zscore=F,
                            gene_weight_method=NULL)
  write.csv(TFA_RNA,paste(result_dir,"TFA_RNA data",date_analysis,".csv",sep = ""))
  
  TFA_RNA_Intron=TFA_RNA[,grep("Intron",colnames(TFA_RNA),value=T)]
  write.csv(TFA_RNA_Intron,paste(result_dir,"TFA_RNA_Intron",date_analysis,".csv",sep = ""))
  
  TFA_RNA_Exon=TFA_RNA[,grep("Exon",colnames(TFA_RNA),value=T)]
  write.csv(TFA_RNA_Exon,paste(result_dir,"TFA_RNA_Exon",date_analysis,".csv",sep = ""))
  
  #calculate circadian related parameters with MetaCycle
  Intron_meta2d_result=meta2d(infile=paste(result_dir,"TFA_RNA_Intron",date_analysis,".csv",sep = ""), filestyle="csv", 
                              outdir=result_dir,
                              timepoints=as.numeric(sub("Intron_","",colnames(TFA_RNA_Intron))),
                              minper = 24,maxper = 24,ARSdefaultPer=24,
                              outputFile=F)
  Intron_meta2d_result=Intron_meta2d_result$meta
  
  Exon_meta2d_result=meta2d(infile=paste(result_dir,"TFA_RNA_Exon",date_analysis,".csv",sep = ""), filestyle="csv", result_dir,
                            timepoints=as.numeric(sub("Exon_","",colnames(TFA_RNA_Exon))),
                            minper = 24,maxper = 24,ARSdefaultPer=24,
                            outputFile=F)
  Exon_meta2d_result=Exon_meta2d_result$meta
  
  Intron_Exon_meta2d_result=merge(Intron_meta2d_result,Exon_meta2d_result,by.x="CycID",by.y="CycID",all=F,
                                  suffixes=c("_Intron","_Exon"))
  write.csv(Intron_Exon_meta2d_result,paste(result_dir,"Intron_Exon_meta2d_result",date_analysis,".csv",sep =""))
  
  #choose TFs with circadian TFA from intron and exon (ARS p<0.05)
  Intron_Exon_meta2d_result_subset=Intron_Exon_meta2d_result[(Intron_Exon_meta2d_result$meta2d_pvalue_Intron<0.05)&
                                                               (Intron_Exon_meta2d_result$meta2d_pvalue_Exon<0.05),]
  Intron_Exon_circadian_TF_vector=Intron_Exon_meta2d_result_subset$CycID
  
  # record the target number of the TFs
  # overlap between dorothea targets and GA_RNA data
  overlap_targets = intersect(rownames(GA_RNA),colnames(dorothea_logic_matrix))
  n_targets = rowSums(dorothea_logic_matrix[,overlap_targets])[Intron_Exon_circadian_TF_vector]
  
  phase_plot_data_frame=Intron_Exon_meta2d_result_subset[,c("CycID","meta2d_phase_Intron","meta2d_phase_Exon"),]
  phase_plot_data_frame = cbind(phase_plot_data_frame,n_targets)
  phase_plot_data_frame$Exon_minus_Intron_phase=phase_plot_data_frame$meta2d_phase_Exon-phase_plot_data_frame$meta2d_phase_Intron 
  
  phase_plot_data_frame$Exon_minus_Intron_phase[phase_plot_data_frame$Exon_minus_Intron_phase<(-12)]=
    phase_plot_data_frame$Exon_minus_Intron_phase[phase_plot_data_frame$Exon_minus_Intron_phase<(-12)]+24
  phase_plot_data_frame$Exon_minus_Intron_phase[phase_plot_data_frame$Exon_minus_Intron_phase>12]=
    phase_plot_data_frame$Exon_minus_Intron_phase[phase_plot_data_frame$Exon_minus_Intron_phase>12]-24
  
  phase_plot_data_frame=phase_plot_data_frame[order(phase_plot_data_frame$Exon_minus_Intron_phase,decreasing = T),]
  phase_plot_data_frame$CycID=factor(phase_plot_data_frame$CycID,levels=phase_plot_data_frame$CycID)
  
}

p1<-ggplot(data=phase_plot_data_frame,mapping=aes(x=CycID,y=Exon_minus_Intron_phase))+
  geom_bar(stat="identity",width=0.5,position='dodge')+
  labs(x = "Exon_minus_Intron_phase", y = "phase difference (h)",
       title = "")+
  theme_bw()+themo_demo+
  geom_text(mapping=aes(x=CycID,y=-1.5,label=paste(n_targets,"\n","targets",sep="")),
            colour="black",position = position_dodge(0),size=5)

pdf(file =paste(result_dir,"phase_difference",date_analysis,".pdf",sep = ""),width = 10,height = 10)
print(p1)
dev.off()

TFA_RNA_Intron_subset=TFA_RNA[Intron_Exon_circadian_TF_vector,grep("Intron",colnames(TFA_RNA),value=T)]
TFA_RNA_Exon_subset=TFA_RNA[Intron_Exon_circadian_TF_vector,grep("Exon",colnames(TFA_RNA),value=T)]

# scale per TF
TFA_RNA_Intron_subset_scale=t(apply(TFA_RNA_Intron_subset,1,scale))
colnames(TFA_RNA_Intron_subset_scale)=colnames(TFA_RNA_Intron_subset)
TFA_RNA_Exon_subset_scale=t(apply(TFA_RNA_Exon_subset,1,scale))
colnames(TFA_RNA_Exon_subset_scale)=colnames(TFA_RNA_Exon_subset)

TFA_all_scale=as.data.frame(cbind(TFA_RNA_Intron_subset_scale,TFA_RNA_Exon_subset_scale))
TFA_all_scale$TF=rownames(TFA_all_scale)
write.csv(TFA_all_scale,paste(result_dir,"TFA_all_scale",".csv",sep = ""))

TFA_all_scale_melt=melt(TFA_all_scale,id.vars = "TF",variable.name = "index",value.name ="scaled_TFA" )
TFA_all_scale_melt$normalization=factor(sub("_.*$","",TFA_all_scale_melt$index),levels=c("Intron","Exon"))
TFA_all_scale_melt$time=as.numeric(sub("^.*_","",TFA_all_scale_melt$index))
TFA_all_scale_melt$TF=factor(TFA_all_scale_melt$TF,levels=phase_plot_data_frame$CycID)
write.csv(TFA_all_scale_melt,paste(result_dir,"TFA_all_scale_melt",".csv",sep = ""))

p_TFA_all_scale_melt<-ggplot(data=TFA_all_scale_melt,mapping=aes(x=time,y=scaled_TFA,colour=normalization))+
  geom_line(size=1)+
  geom_point(size = 2)+
  labs(x = "Time (h)", y = "scaled_TFA",
       title ="")+
  scale_x_continuous(limits=c(0,22),breaks=seq(0,22,2))+
  facet_grid(TF~.)+
  theme_bw()+themo_demo

pdf(paste(result_dir,"p_TFA_all_scale_melt",date_analysis,".pdf",sep = ""),height=12,width=5.5)
print(p_TFA_all_scale_melt)
dev.off()


# analysis for expression data
# change MetaCycle input: TF activity -> TF expression (RPKM matrix)
GA_RNA_exon = GA_RNA[,13:24]
match_id = match(Intron_Exon_circadian_TF_vector,rownames(GA_RNA_exon))
names(match_id) = Intron_Exon_circadian_TF_vector
# not all TFs have expression level in the matrix
GA_RNA_exon_sub = GA_RNA_exon[match_id,]
rownames(GA_RNA_exon_sub) = Intron_Exon_circadian_TF_vector
# scale
TF_exon_scale =t(apply(GA_RNA_exon_sub,1,scale))
colnames(TF_exon_scale)=colnames(GA_RNA_exon_sub)

TF_exon_scale_melt=melt(TF_exon_scale,id.vars = "TF",variable.name = "index",value.name ="scaled_TF_exon" )
colnames(TF_exon_scale_melt) = c("TF","index","scaled_TF_exon")
TF_exon_scale_melt$normalization=factor(sub("_.*$","",TF_exon_scale_melt$index),levels=c("Intron","Exon"))
TF_exon_scale_melt$time=as.numeric(sub("^.*_","",TF_exon_scale_melt$index))
TF_exon_scale_melt$TF=factor(TF_exon_scale_melt$TF,levels=phase_plot_data_frame$CycID)
TF_exon_scale_melt$normalization = "Exon_exp"

# plot dynamics of TF expression
p_TF_exon_scale_melt<-ggplot(data=TF_exon_scale_melt,mapping=aes(x=time,y=scaled_TF_exon,colour=normalization))+
  geom_line(size=1)+
  geom_point(size = 2)+
  labs(x = "Time (h)", y = "scaled_TF_exon",
       title ="")+
  scale_x_continuous(limits=c(0,22),breaks=seq(0,22,2))+
  facet_grid(TF~.)+
  theme_bw()+themo_demo

pdf(paste(result_dir,"p_TF_exon_scale_melt",date_analysis,".pdf",sep = ""),height=12,width=5.5)
print(p_TF_exon_scale_melt)
dev.off()

# Compare TF activity with TF expression
TF_exon_scale_melt$normalization = "Exon_exp"
colnames(TF_exon_scale_melt) = colnames(TFA_all_scale_melt)
TFA_all_scale_melt_both = rbind(TF_exon_scale_melt,TFA_all_scale_melt)
TFA_all_scale_melt_both$normalization=factor(TFA_all_scale_melt_both$normalization,levels=c("Intron","Exon","Exon_exp"))

p_TFA_all_scale_melt_both<-ggplot(data=TFA_all_scale_melt_both,mapping=aes(x=time,y=scaled_TFA,colour=normalization))+
  geom_line(size=1)+
  geom_point(size = 2)+
  labs(x = "Time (h)", y = "scaled_TFA",
       title ="")+
  scale_x_continuous(limits=c(0,22),breaks=seq(0,22,2))+
  facet_grid(TF~.)+
  theme_bw()+themo_demo

pdf(paste(result_dir,"p_TFA_all_scale_melt_both",date_analysis,".pdf",sep = ""),height=12,width=5.5)
print(p_TFA_all_scale_melt_both)
dev.off()

# calculate correlation between TFA and expression
dim(TFA_RNA_Exon)
dim(TFA_RNA_Intron)
dim(GA_RNA_exon)
genes.overlap = intersect(rownames(GA_RNA_exon),rownames(TFA_RNA_Exon))
TF_exp_over = GA_RNA_exon[genes.overlap,]
TFA_in_over = TFA_RNA_Intron[genes.overlap,]
TFA_ex_over = TFA_RNA_Exon[genes.overlap,]
# dataframe of TFA-expression correlation
df.cor = NULL
for(i in 1:nrow(TF_exp_over)){
  cor.in = cor(TF_exp_over[i,],TFA_in_over[i,])
  cor.ex = cor(TF_exp_over[i,],TFA_ex_over[i,])
  df.cor = rbind(df.cor,c(cor.in,cor.ex))
}
rownames(df.cor) = rownames(TF_exp_over)
colnames(df.cor) = c("cor.in_exp","cor.ex_exp")
df.cor = as.data.frame(df.cor)
hist(df.cor$cor.in_exp)
hist(df.cor$cor.ex_exp)
plot(df.cor$cor.ex_exp,df.cor$cor.in_exp)
abline(a = 0, b = 1, col = "green")
plot(abs(df.cor$cor.ex_exp),abs(df.cor$cor.in_exp))
abline(a = 0, b = 1, col = "green")
# select circadian TFs
df.circadian = data.frame(meta2d_pvalue_Exon = Intron_Exon_meta2d_result$meta2d_pvalue_Exon,
                          meta2d_pvalue_Intron = Intron_Exon_meta2d_result$meta2d_pvalue_Intron)
rownames(df.circadian) = Intron_Exon_meta2d_result$CycID
hist(df.circadian$meta2d_pvalue_Exon)
# dataframe of all information
df.summary = merge(df.circadian,df.cor,by = "row.names",all= T)
# circadian TFs
pdf(paste(result_dir,"Select Circadian TFs",".pdf",sep = ""),height=9,width=9)

plot(-log10(df.summary$meta2d_pvalue_Exon),-log10(df.summary$meta2d_pvalue_Intron),main = "Select Circadian TFs",
     xlab = "-log(p-value) exon TFA",ylab = "-log(p-value) intron TFA")
abline(h = -log10(0.05),col = "red",lty = 2)
abline(v = -log10(0.05),col = "blue",lty = 2)
thr = 0.05
nn = which(df.summary$meta2d_pvalue_Exon>thr&df.summary$meta2d_pvalue_Intron>thr)
pn = which(df.summary$meta2d_pvalue_Exon<thr&df.summary$meta2d_pvalue_Intron>thr)
np = which(df.summary$meta2d_pvalue_Exon>thr&df.summary$meta2d_pvalue_Intron<thr)
pp = which(df.summary$meta2d_pvalue_Exon<thr&df.summary$meta2d_pvalue_Intron<thr)
points(-log10(df.summary$meta2d_pvalue_Exon[nn]),-log10(df.summary$meta2d_pvalue_Intron[nn]),col = "black",pch = 16)
points(-log10(df.summary$meta2d_pvalue_Exon[pn]),-log10(df.summary$meta2d_pvalue_Intron[pn]),col = "blue",pch = 16)
points(-log10(df.summary$meta2d_pvalue_Exon[np]),-log10(df.summary$meta2d_pvalue_Intron[np]),col = "red",pch = 16)
points(-log10(df.summary$meta2d_pvalue_Exon[pp]),-log10(df.summary$meta2d_pvalue_Intron[pp]),col = "purple",pch = 16)
pie(c(length(nn),length(pn),length(pp),length(np)),col = c("black","blue","purple","red"),
    labels = c(length(nn),length(pn),length(pp),length(np)))
dev.off()

circadian_TFs = as.character(df.summary$Row.names[pp])
# show all genes
nn
df.summary$Row.names[nn]
df.summary$Row.names[pn]
df.summary$Row.names[np]
df.summary$Row.names[pp]
# compare circadian genes with non-circadian genes
if(T){
  # load circadian genes
  # GO: circadian clock
  GO_CC = read.table(file = "../../materials/circadian_clock_genes_GO.txt",sep = "\t",header = F)
  GO_CC = GO_CC$V2
  intersect(GO_CC,df.summary$Row.names)
  library(VennDiagram)
  venn.diagram(x = list(GO_CC,df.summary$Row.names),category.names = c("Circadian Clock Genes (GO)","Dorothea TFs"),
               filename = "../results/venn.png",height = 900, width = 1600,
               cat.default.pos = "text",cat.cex = 0.5,cat.pos = c(0,0),cat.dist = c(0.2, 0.3))
  # save gene names
  # write.table(df.summary$Row.names,file = "genenames.txt",row.names = F,col.names = F,quote = F,sep = "\t")
  # circadian genes from 2014,Cell, BinFang et al, S2
  geneset = c("BMAL1","CLOCK","NPAS2","FOXA1","HNF4A","HDAC3")
  # circadian genes from 2016,literature
  
  genes = read.table("../results/gene_names.txt")
  art_num = read.table("../results/article_numbers.txt")
  df.art_num = data.frame(gene = as.character(genes), art_num = as.numeric(art_num))
  plot(1:nrow(df.art_num),sort(df.art_num$art_num,decreasing = T),log = "y", ylab = "number of articles", xlab = "Dorothea TFs",
       main = "Relation to Circadian Clock")
  # df.art_num = df.art_num[order(df.art_num$art_num,decreasing = T),]
  
  df.summary$art_num = df.art_num$art_num[match(df.summary$Row.names,df.art_num$gene)]
  plot(df.summary$meta2d_pvalue_Exon,df.summary$meta2d_pvalue_Intron,log = "xy",main = "Select Circadian TFs")
  abline(h = 0.05,col = "red")
  abline(v = 0.05,col = "blue")
  id = which(df.summary$art_num>=5)
  points(df.summary$meta2d_pvalue_Exon[id],df.summary$meta2d_pvalue_Intron[id],col = "green",pch = 16)
  
  length(intersect(nn,id))/length(nn)
  length(intersect(pn,id))/length(pn)
  length(intersect(np,id))/length(np)
  length(intersect(pp,id))/length(pp)
  
  df.summary$Row.names[intersect(np,id)]
  
  # plot threshold - overlap curve
  df.thr_over = NULL
  for(thr in sort(unique(df.summary$art_num))){
    id = which(df.summary$art_num>=thr)
    df.thr_over = rbind(df.thr_over,c(thr,length(id),length(intersect(nn,id))/length(nn),length(intersect(pn,id))/length(pn),
                                      length(intersect(np,id))/length(np),length(intersect(pp,id))/length(pp),
                                      length(intersect(unique(c(pn,pp)),id))/length(unique(c(pn,pp))),
                                      length(intersect(unique(c(np,pp)),id))/length(unique(c(np,pp)))))
  }
  colnames(df.thr_over) = c("thr","num_cc","p.nn","p.pn","p.np","p.pp","p.ex","p.in")
  df.thr_over  = as.data.frame(df.thr_over)
  plot(df.thr_over$thr,df.thr_over$p.nn,log = "x",type = "l",col = "black",xlab = "article number",ylab = "Circadian proportion",
       main = "Compare four genesets")
  points(df.thr_over$thr,df.thr_over$p.pn,type = "l",col = "blue")
  points(df.thr_over$thr,df.thr_over$p.np,type = "l",col = "red")
  points(df.thr_over$thr,df.thr_over$p.pp,type = "l",col = "purple")
  plot(df.thr_over$thr,df.thr_over$p.nn,log = "x",type = "l",col = "black",xlab = "article number",ylab = "Circadian proportion",
       main = "Compare intron and exon TFA")
  points(df.thr_over$thr,df.thr_over$p.ex,type = "l",col = "blue")
  points(df.thr_over$thr,df.thr_over$p.in,type = "l",col = "red")
}


# plot TFA - expression correlation
plot(df.summary$cor.ex_exp, df.summary$meta2d_pvalue_Exon,log = "y",main = "TFA - expression correlation (exon)")
points(df.summary$cor.ex_exp[which(df.summary$meta2d_pvalue_Exon>0.05)],
       df.summary$meta2d_pvalue_Exon[which(df.summary$meta2d_pvalue_Exon>0.05)],col ="grey")
abline(v = 0.5,lty = 2)
abline(v = -0.5,lty = 2)
abline(h = 0.05,col = "grey",lty = 2)

plot(df.summary$cor.in_exp, df.summary$meta2d_pvalue_Intron,log = "y",main = "TFA - expression correlation (intron)")
points(df.summary$cor.in_exp[which(df.summary$meta2d_pvalue_Intron>0.05)],
       df.summary$meta2d_pvalue_Intron[which(df.summary$meta2d_pvalue_Intron>0.05)],col ="grey")
abline(v = 0.5,lty = 2)
abline(v = -0.5,lty = 2)
abline(h = 0.05,col = "grey",lty = 2)

# intron - exon correlation
plot(df.summary$cor.ex_exp,df.summary$cor.in_exp,main = "cor(TFA,exp): intron vs exon")
points(df.summary$cor.ex_exp[pp],df.summary$cor.in_exp[pp],col = "purple",pch = 16)
points(df.summary$cor.ex_exp[pn],df.summary$cor.in_exp[pn],col = "blue",pch = 16)
points(df.summary$cor.ex_exp[np],df.summary$cor.in_exp[np],col = "red",pch = 16)
abline(v = 0.5,lty = 2)
abline(v = -0.5,lty = 2)
abline(h = 0.5,lty = 2)
abline(h = -0.5,lty = 2)

# investigate circadian TFs
plot(df.summary$cor.ex_exp[pp],df.summary$cor.in_exp[pp],col = "purple",pch = 16,main = "cor(TFA,exp): intron vs exon (circadian TFs)")
abline(v = 0.5,lty = 2)
abline(v = -0.5,lty = 2)
abline(h = 0.5,lty = 2)
abline(h = -0.5,lty = 2)
# correlation for circadian TFs
df.cor.cir = df.cor[circadian_TFs,]
rownames(df.cor.cir) = circadian_TFs
# intron
barplot(sort(t(df.cor.cir)[1,]),ylim = c(-1,1),main = "TFA(intron) - expression correlation")
abline(h = 0.5,lty = 2)
abline(h = -0.5,lty = 2)
# exon
barplot(sort(t(df.cor.cir)[2,]),ylim = c(-1,1),main = "TFA(exon) - expression correlation")
abline(h = 0.5,lty = 2)
abline(h = -0.5,lty = 2)

# ggplot
library(ggplot2)
df.plot = data.frame(TF = rep(rownames(df.cor.cir),2), cor = c(df.cor.cir$cor.in_exp,df.cor.cir$cor.ex_exp),
                     data_type = c(rep("intron",nrow(df.cor.cir)),rep("exon",nrow(df.cor.cir))))
df.plot = df.plot[order(df.plot$cor),]
df.plot$TF = factor(df.plot$TF,levels = unique(df.plot$TF[which(df.plot$data_type=="intron")]))
df.plot$data_type = factor(df.plot$data_type,levels = c("intron","exon"))
p = ggplot() + geom_bar(data = df.plot, aes(x = TF,y = cor, fill = data_type),stat = "identity",
                        position = "dodge")+ scale_fill_manual(values=c("#F8766D","#00BFC4","grey"))

p = p+ labs(x = NULL,y = "correlation",title = "TFA - expression correlation")
p = p+ geom_hline(aes(yintercept = 0.5),linetype= "dashed") + geom_hline(aes(yintercept = -0.5),linetype= "dashed") + ylim(-1,1)
p

# plot correlation difference - phase difference
cor.diff = df.summary$cor.in_exp[pp] - df.summary$cor.ex_exp[pp]
names(cor.diff) = df.summary$Row.names[pp]
phase.diff = phase_plot_data_frame$Exon_minus_Intron_phase
names(phase.diff) = phase_plot_data_frame$CycID
cor.diff = cor.diff[match(names(phase.diff),names(cor.diff))]
plot(phase.diff,cor.diff,main = "Phase difference vs Correlation difference")
abline(h = 0,lty = 2)

# plot phase difference - p-value

# investigate the intron length of targets
if(T){
  # intron length information
  load(file = "../../materials/Intron_length_EnsDb.Mmusculus.v79.RData")
  df.intron_length  = df
  
  # extract TF - targets
  targets = dorothea_logic_matrix[Intron_Exon_circadian_TF_vector,overlap_targets]
  # transfer matrix to linklist
  mat_to_linklist = function(a){
    TFs = rownames(a)
    targets = colnames(a)
    linklist = NULL
    for(i in 1:nrow(a)){
      for(j in 1:ncol(a)){
        if(a[i,j]==1){
          linklist = rbind(linklist,c(TFs[i],targets[j]))
        }
      }
    }
    linklist = as.data.frame(linklist)
    colnames(linklist) = c("TF","target")
    return(linklist)
  }
  linklist = mat_to_linklist(targets)
  
  # intron length information and time delay information
  match(linklist$target,df.intron_length$gene)
  phase_plot_data_frame
  # transfer TF-target linklist to timedelay - intron length linklist
  linklist.time_intronlen = data.frame(timedelay = phase_plot_data_frame$Exon_minus_Intron_phase[match(linklist$TF,phase_plot_data_frame$CycID)],
                                       intron_length = df.intron_length$intron_length[match(linklist$target,df.intron_length$gene)])
  x = linklist.time_intronlen$timedelay
  y = log10(linklist.time_intronlen$intron_length+1)
  pdf(paste(result_dir,"Circadian Effect of Intron Length",".pdf",sep = ""),height=9,width=9)
  plot(x,y,main = "Intron Length - Phase difference",pch = 16,
       xlab = "Phase difference of TF (h)", ylab = "log10(intron length of targets + 1)")
  dev.off()
  cor(x,y)
}

# investigate the mRNA half-life of targets

if(T){
  # load mRNA half-life data
  if(T){
    hl = read.csv(file = "../../materials/2009_mESC.csv")
    non_blank_id =which(hl$Gene.symbol != '')
    hl = hl[non_blank_id,]
    # remove zero half-life genes
    non_zero_genes = which(hl$Halflife.for.all.cells!=0)
    hl = hl[non_zero_genes,]
    
    hl.data = hl$Halflife.for.all.cells
    names(hl.data) = hl$Gene.symbol
    hl = hl.data
  }
  
  # extract TF - targets
  targets = dorothea_logic_matrix[Intron_Exon_circadian_TF_vector,overlap_targets]
  # transfer matrix to linklist
  mat_to_linklist = function(a){
    TFs = rownames(a)
    targets = colnames(a)
    linklist = NULL
    for(i in 1:nrow(a)){
      for(j in 1:ncol(a)){
        if(a[i,j]==1){
          linklist = rbind(linklist,c(TFs[i],targets[j]))
        }
      }
    }
    linklist = as.data.frame(linklist)
    colnames(linklist) = c("TF","target")
    return(linklist)
  }
  linklist = mat_to_linklist(targets)
  
  # mRNA half-life information and time delay information
  match(linklist$target,names(hl))
  phase_plot_data_frame
  # transfer TF-target linklist to timedelay - intron length linklist
  linklist.time_mRNAhl = data.frame(timedelay = phase_plot_data_frame$Exon_minus_Intron_phase[match(linklist$TF,phase_plot_data_frame$CycID)],
                                       mRNA_hl = hl[match(linklist$target,names(hl))])
  nrow(na.omit(linklist.time_mRNAhl))
  nrow((linklist.time_mRNAhl))
  linklist.time_mRNAhl = na.omit(linklist.time_mRNAhl)
  x = linklist.time_mRNAhl$timedelay
  y = linklist.time_mRNAhl$mRNA_hl
  plot(x,y,main = "mRNA half-life - time delay",
       xlab = "time delay of TF (h)", ylab = "mRNA half-life (h)")
  cor(x,y)
}

# calculate TFA for additional TFs
if(T){
  # linklist to matrix function
  link_to_mat = function(a){
    TFs = unique(as.character(as.matrix(a[,1])))
    targets = unique(as.character(as.matrix(a[,3])))
    mat = matrix(0,nrow = length(TFs),ncol = length(targets))
    rownames(mat) = TFs
    colnames(mat) = targets
    for(i in 1:nrow(a)){
      mat[as.character(a$tf[i]),as.character(a$target[i])] = 1
    }
    return(mat)
  }
  load(file ="C:/Users/HP/Documents/dorothea/dorothea_database.RData")
  TFs = c("Clock","Arntl","Npas2","Nfil3","Nr1d1","Nr1d2","Rora","Rorb")
  ids = which(is.element(dorothea_mm$tf,TFs)&dorothea_mm$mor==1)
  sublinks = dorothea_mm[ids,]
  # TFA for each TF
  dorothea_sub = link_to_mat(sublinks)
  TFA_sub=scTFA_calculation(as.matrix(GA_RNA),dorothea_sub,zscore=F,gene_weight_method=NULL)
  
  par(mfrow =c(2,2))
  for(i in 1:nrow(TFA_sub)){
    plot(seq(0,22,2),scale(TFA_sub[i,1:12]),col = "red",type = "b",ylim = c(-2,2),main = rownames(TFA_sub)[i],
         xlab = "time", ylab = "z-score")
    points(seq(0,22,2),scale(TFA_sub[i,13:24]),col = "blue",type = "b")
  }
  
  par(mfrow =c(2,2))
  for(i in 1:nrow(TFA_sub)){
    plot(seq(0,22,2),scale(TFA_sub[i,1:12]),col = "red",type = "b",ylim = c(-2,2),main = rownames(TFA_sub)[i],
         xlab = "time", ylab = "z-score")
    points(seq(0,22,2),scale(TFA_sub[i,13:24]),col = "blue",type = "b")
  }
  
  pdf(paste(result_dir,"Circadian TFA dynamics",".pdf",sep = ""),height=9,width=9)
  par(mfrow =c(1,1))
  i=1
  plot(seq(0,22,2),scale(TFA_sub[i,1:12]),col = "red",type = "b",ylim = c(-2,2),main = paste(rownames(TFA_sub)[i],"Dynamics"),
       xlab = "Time (h)", ylab = "TFA Level")
  points(seq(0,22,2),scale(TFA_sub[i,13:24]),col = "blue",type = "b")
  i=2
  plot(seq(0,22,2),scale(TFA_sub[i,1:12]),col = "red",type = "b",ylim = c(-2,2),main =  paste(rownames(TFA_sub)[i],"Dynamics"),
       xlab = "Time (h)", ylab = "TFA Level")
  points(seq(0,22,2),scale(TFA_sub[i,13:24]),col = "blue",type = "b")
  dev.off()
  
  
  
}

# plot for phase difference - robustness
pdf(paste(result_dir,"Robustness of results",".pdf",sep = ""),height=9,width=9)
if(T){
  phase_plot_data_frame$Freq = robustness$Freq[match(phase_plot_data_frame$CycID,robustness$Var1)]
  plot(phase_plot_data_frame$Exon_minus_Intron_phase,phase_plot_data_frame$Freq/4,ylim = c(0,1),xlab = "phase difference (h)",
       ylab = "Robustness",main = "Show robustness of results")
  abline(v = 0,col = "black",lty = 2)
  abline(h = 0.5,col = "red",lty = 2)
}
dev.off()
