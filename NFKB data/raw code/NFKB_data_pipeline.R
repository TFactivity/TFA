
library(Seurat)
require(ggplot2)
require(reshape2)
library("org.Mm.eg.db")
library(EnsDb.Mmusculus.v79)
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
  axis.text.y = element_text(size=22)
)

data_dir="../raw data/"
result_dir="../results/"

#import RNA data
RNA_fpkm=read.csv(paste(data_dir,"GSE94383_fpms.csv",sep=""),row.names = 1)
colnames(RNA_fpkm)=sub("\\.","-",colnames(RNA_fpkm))

#import metadata
RNA_cell_id=read.csv(paste(data_dir,"GSE94383_cell_ids.csv",sep=""),row.names = 1)
RNA_cell_id=RNA_cell_id[RNA_cell_id$Condition %in% c("NoStim","Stim"),]
rownames(RNA_cell_id)=paste("X",rownames(RNA_cell_id),sep="")

#import movie data
RNA_movie_data=read.csv(paste(data_dir,"GSE94383_single_cell_dynamics.csv",sep=""),row.names = 1)
rownames(RNA_movie_data)=paste("X",rownames(RNA_movie_data),sep="")
RNA_movie_data$movie_time_points=apply(RNA_movie_data,1,function(x) sum(!is.na(x)))

#create seurat object
RAW_seurat=CreateSeuratObject(counts = RNA_fpkm, project = "RAW_fpkm", min.cells = 0, min.features = 0)

#add info of movie and cell id to metadata
#View(A549_seurat@meta.data)
RAW_seurat@meta.data=cbind(RAW_seurat@meta.data,RNA_cell_id[rownames(RAW_seurat@meta.data),],
                           RNA_movie_data[rownames(RAW_seurat@meta.data),])

#import GRN data
mouse_dorothea_logic_matrix=read.csv("../../materials/mouse_GRN_classAB_activation_only/dorothea_mouse_AB.csv",
                                     row.names = 1, check.names = F)

#Convert from gene.symbol to ensembl.gene
geneSymbols <-colnames(mouse_dorothea_logic_matrix)
geneIDs <- ensembldb::select(EnsDb.Mmusculus.v79, keys= geneSymbols, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))

mouse_dorothea_logic_matrix_ensembl=matrix(0,nrow(mouse_dorothea_logic_matrix),1,
                                           dimnames = list(rownames(mouse_dorothea_logic_matrix),"NA"))
for(col_i in colnames(mouse_dorothea_logic_matrix)){
  if(col_i %in% geneIDs$SYMBOL){
    gene_id_i=geneIDs[geneIDs$SYMBOL==col_i,"GENEID"]
    for(j_gene_id_i in gene_id_i){
      mouse_dorothea_logic_matrix_ensembl=cbind(mouse_dorothea_logic_matrix_ensembl,mouse_dorothea_logic_matrix[,col_i])
      colnames(mouse_dorothea_logic_matrix_ensembl)[ncol(mouse_dorothea_logic_matrix_ensembl)]=j_gene_id_i
    }
  }
}
mouse_dorothea_logic_matrix_ensembl=mouse_dorothea_logic_matrix_ensembl[,-1]


#import intron/exon CPM data and compute TFA
out_file_dir="../raw data/tpmcalculator_out_file/"
out_file_name_vector=list.files(out_file_dir)
out_file_data=list()

# 20 min
i = 1
for(out_file_i in out_file_name_vector){
  data_i=read.table(paste(out_file_dir,out_file_i,sep=""),fill = T,header=T)
  data_i=na.omit(data_i)
  rownames(data_i)=sub("\\..*$","",data_i$Gene_Id)
  
  data_i$totalCPM=1e6*(data_i$Reads)/sum(data_i$Reads)
  data_i$ExonCPM=1e6*(data_i$ExonReads)/sum(data_i$Reads)
  data_i$IntronCPM=1e6*(data_i$IntronReads)/sum(data_i$Reads)
  
  TFA_ExonCPM=as.data.frame(scTFA_calculation(as.matrix(data_i[,"ExonCPM",drop=F]),mouse_dorothea_logic_matrix_ensembl,zscore=F,
                                                                 gene_weight_method=NULL))
  colnames(TFA_ExonCPM)="TFA_ExonCPM"
  TFA_IntronCPM=as.data.frame(scTFA_calculation(as.matrix(data_i[data_i$IntronLength>0,"IntronCPM",drop=F]),mouse_dorothea_logic_matrix_ensembl,zscore=F,
                                                                   gene_weight_method=NULL))
  colnames(TFA_IntronCPM)="TFA_IntronCPM"
  
  TFA_all=cbind(TFA_IntronCPM[rownames(TFA_IntronCPM),,drop=F],
                TFA_ExonCPM[rownames(TFA_IntronCPM),,drop=F])
  TFA_all$TF_name=rownames(TFA_all)
  out_file_data[[out_file_i]]=TFA_all
  
  print(i/length(out_file_name_vector))
  i = i + 1
}

#change names
cell_sra_result=read.csv("../raw data/sra_result.csv")
cell_SraRunInfo=read.csv("../raw data/SraRunInfo.csv")
cell_metadata=merge(cell_sra_result,cell_SraRunInfo,all=F,by.x="Experiment.Accession",by.y="Experiment")
rownames(cell_metadata)=cell_metadata$Run

names_out_file_data_new=sub("_genes.out","",names(out_file_data))
names_out_file_data_new=cell_metadata[names_out_file_data_new,"Experiment.Title"]
names_out_file_data_new=sub("\\).*$","",sub("^.*\\(","",names_out_file_data_new))
names_out_file_data_new=paste("X",names_out_file_data_new,sep="")

names(out_file_data)=names_out_file_data_new
#saveRDS(out_file_data,file=paste(result_dir,"out_file_data",".rds",sep = ""))
out_file_data=readRDS(paste(result_dir,"out_file_data.rds",sep=""))

#subset
time_point_vector=c(75,150,300)
max_movie_time_point_vector=c(75,150,300)/5+1
TF_chosen=c("Rela")
# number of targets
rowSums(mouse_dorothea_logic_matrix)["Rela"]
TF_targets = colnames(mouse_dorothea_logic_matrix)[which(mouse_dorothea_logic_matrix["Rela",]==1)]
#overlap with RNA seq data

correlation_matrix=matrix(,length(time_point_vector),3,
                          dimnames = list(paste("time",time_point_vector,sep="_"),c("IntronCPM_vs_movie","ExonCPM_vs_movie","cellnum")))
for(i_time_point in 1:length(time_point_vector)){
  time_point=time_point_vector[i_time_point]
  max_movie_time_point=max_movie_time_point_vector[i_time_point]
  RAW_subset=subset(RAW_seurat, subset = Time.point == time_point)
  
  RAW_subset=subset(RAW_subset, subset = nCount_RNA < 9e6 ) ##according to ???###
  
  out_file_data_subset=out_file_data[which(names(out_file_data) %in% colnames(RAW_subset))]
  RAW_subset=RAW_subset[,intersect(colnames(RAW_subset),names(out_file_data_subset))]
  
  #merge data
  index=0
  for(name_i in names(out_file_data_subset)){
    colnames(out_file_data_subset[[name_i]])=paste(colnames(out_file_data_subset[[name_i]]),name_i,sep="_")
    if(index==0){
      out_file_data_merged=out_file_data_subset[[name_i]]
    }else{
      out_file_data_merged=merge(out_file_data_merged,out_file_data_subset[[name_i]],all=F,
                                 by.x=paste("TF_name",names(out_file_data_subset)[1],sep="_"),
                                 by.y=paste("TF_name",name_i,sep="_"))
    }
    index=index+1
  }
  rownames(out_file_data_merged)=out_file_data_merged[,paste("TF_name",names(out_file_data_subset)[1],sep="_")]
  write.csv(out_file_data_merged,paste("out_file_data_merged",time_point,".csv"))
  
  TFA_ExonCPM_all=out_file_data_merged[,paste("TFA_ExonCPM",colnames(RAW_subset),sep="_")]
  colnames(TFA_ExonCPM_all)=colnames(RAW_subset)
  rownames(TFA_ExonCPM_all)=paste("TFA_ExonCPM",rownames(TFA_ExonCPM_all),sep="_")
  
  TFA_IntronCPM_all=out_file_data_merged[,paste("TFA_IntronCPM",colnames(RAW_subset),sep="_")]
  colnames(TFA_IntronCPM_all)=colnames(RAW_subset)
  rownames(TFA_IntronCPM_all)=paste("TFA_IntronCPM",rownames(TFA_IntronCPM_all),sep="_")
  
  #add metadata
  movie_TFA_data=cbind(RAW_subset@meta.data,
                       t(TFA_ExonCPM_all)[rownames(RAW_subset@meta.data),],
                       t(TFA_IntronCPM_all)[rownames(RAW_subset@meta.data),])
  
  index_vector=c("TFA_IntronCPM","TFA_ExonCPM")
  
  #correlation between TFA and nucl loc of every frame
  row_names_every=paste("corr with nucl of the",1:(max_movie_time_point+1),"frame before")
  corr_matrix_every=matrix(NA,length(row_names_every),length(index_vector),
                           dimnames = list(row_names_every,index_vector))
  for(index_i in index_vector){
    for(time_point_num in 1:(max_movie_time_point+1)){
      corr_matrix_every[paste("corr with nucl of the",time_point_num,"frame before"),index_i]=
        cor(movie_TFA_data[,paste(index_i,TF_chosen,sep="_")],
            movie_TFA_data[,paste("X",max_movie_time_point-time_point_num+1,sep="")],use="pairwise.complete.obs")
    }
  }
  corr_matrix_every=as.data.frame(corr_matrix_every)
  corr_matrix_every$label=as.numeric(sub(" frame before","",sub("corr with nucl of the ","",rownames(corr_matrix_every))))
  # write.csv(corr_matrix_every,paste("corr_matrix_every data",time_point,".csv"))
  
  correlation_matrix[paste("time",time_point,sep="_"),c("IntronCPM_vs_movie","ExonCPM_vs_movie","cellnum")]=
    c(as.numeric(corr_matrix_every[corr_matrix_every$label==1,c("TFA_IntronCPM","TFA_ExonCPM")]),ncol(RAW_subset))
  
  #correlation between TFA and nucl loc of last frames
  row_names_last=paste("corr with sum nucl last",1:(max_movie_time_point+1),"frame",sep="_")
  corr_matrix_last=matrix(NA,length(row_names_last),length(index_vector),
                          dimnames = list(row_names_last,index_vector))
  for(index_i in index_vector){
    for(time_point_num in 1:(max_movie_time_point+1)){
      corr_matrix_last[paste("corr with sum nucl last",time_point_num,"frame",sep="_"),index_i]=
        cor(movie_TFA_data[,paste(index_i,TF_chosen,sep="_")],
            apply(movie_TFA_data[,paste("X",(max_movie_time_point-time_point_num+1):max_movie_time_point,sep=""),drop=F],1,sum),
            use="pairwise.complete.obs")
    }
  }
  corr_matrix_last=as.data.frame(corr_matrix_last)
  corr_matrix_last$label=sub("_frame","",sub("corr with sum nucl last_","",rownames(corr_matrix_last)))
  corr_matrix_last$label=as.numeric(corr_matrix_last$label)
  # write.csv(corr_matrix_last,paste("corr_matrix_last data",time_point,".csv"))
}
correlation_matrix=as.data.frame(correlation_matrix)
correlation_matrix$time=factor(rownames(correlation_matrix),levels=rownames(correlation_matrix))
# write.csv(correlation_matrix,paste("correlation_matrix data",".csv"))

correlation_matrix_plot=melt(correlation_matrix,id.vars = c("time","cellnum"),variable.name = "index",value.name = "corr")
# write.csv(correlation_matrix_plot,paste("correlation_matrix_plot",".csv"))

p1<-ggplot(data=correlation_matrix_plot,mapping=aes(x=time,y=corr,fill=index,group=index))+
  geom_bar(stat="identity",position='dodge',width=0.5)+
  labs(x = "time", y = "Pearson correlation",
       title = paste("TF",TF_chosen))+
  geom_text(mapping=aes(x=time,y=corr+0.05,label=round(corr,3)),colour="black",position = position_dodge(0.5),size=6)+
  geom_text(mapping=aes(x=time,y=0,label=paste(cellnum,"\n","cell",sep="")),colour="black",position = position_dodge(0.5),size=6)+
  scale_y_continuous(limits=c(-0.03,0.53),breaks = seq(0,0.5,0.1))+
  theme_bw()+themo_demo

pdf(file =paste(result_dir,"correlation_matrix_plot",TF_chosen,".pdf",sep = ""),width = 10,height = 8)
print(p1)
dev.off()

# compare with TF expression
# TFA data for each cell
TF = "Rela"
TFA_cells = NULL
for(i in 1:length(out_file_data)){
  TFA_cells = rbind(TFA_cells,out_file_data[[i]][TF,1:2])
}
rownames(TFA_cells) = names(out_file_data)
TFE_cells = NULL
# load TF expression data
if(F){
  out_file_dir="../raw data/tpmcalculator_out_file/"
  out_file_name_vector=list.files(out_file_dir)
  out_file_data=list()
  
  # 20 min
  i = 1
  for(out_file_i in out_file_name_vector){
    data_i=read.table(paste(out_file_dir,out_file_i,sep=""),fill = T,header=T)
    data_i=na.omit(data_i)
    rownames(data_i)=sub("\\..*$","",data_i$Gene_Id)
    
    data_i$totalCPM=1e6*(data_i$Reads)/sum(data_i$Reads)
    data_i$ExonCPM=1e6*(data_i$ExonReads)/sum(data_i$Reads)
    data_i$IntronCPM=1e6*(data_i$IntronReads)/sum(data_i$Reads)
    
    geneID <- ensembldb::select(EnsDb.Mmusculus.v79, keys= TF, keytype = "SYMBOL", columns = c("SYMBOL","GENEID"))
    TFE_cells = c(TFE_cells,data_i[geneID$GENEID,"ExonCPM"])
    
    print(i/length(out_file_name_vector))
    i = i + 1
  }
  
}
# save(TFE_cells,file = paste(result_dir,"TFE_cells.RData",sep = ""))
load(paste(result_dir,"TFE_cells.RData",sep = ""))
df.TFA = data.frame(TFA.in = TFA_cells$TFA_IntronCPM, TFA.ex = TFA_cells$TFA_ExonCPM, TFE = TFE_cells)
rownames(df.TFA) = rownames(TFA_cells)
cor(na.omit(df.TFA))
nrow(na.omit(df.TFA))
nrow(df.TFA)
# calculate total TFA
df.TFA$TFA.to = df.TFA$TFA.in+df.TFA$TFA.ex
# metadata
tp = c(0,75,150,300)
df.all = cbind(log2(df.TFA+1),RNA_cell_id$Time.point[na.omit(match(rownames(df.TFA),rownames(RNA_cell_id)))])
colnames(df.all)[5] = "time"
df.all$time = as.numeric(as.character(df.all$time))
df.all$time = factor(df.all$time,levels = c(0,75,150,300))
table(df.all$time)
# plot(df.all$time,df.all$TFE,log ="y")

# plot
library(ggplot2)
p1 = ggplot(df.all,aes(x = time,y = TFE)) + geom_violin(aes(fill = time)) + ylim(0, 10)
p2 = ggplot(df.all,aes(x = time,y = TFA.in))+geom_violin(aes(fill = time)) + ylim(0, 10)
p3 = ggplot(df.all,aes(x = time,y = TFA.ex))+geom_violin(aes(fill = time)) + ylim(0, 10)
library("easyGgplot2")
ggplot2.multiplot(p1,p2,p3,cols = 3)

# Rela activity data
time_point_vector=c(75,150,300)
max_movie_time_point_vector=time_point_vector/5+1
i_time_point = 2
for(i_time_point in 1:length(time_point_vector)){
  time_point=time_point_vector[i_time_point]
  max_movie_time_point=max_movie_time_point_vector[i_time_point]
  RAW_subset=subset(RAW_seurat, subset = Time.point == time_point)
  RAW_subset=subset(RAW_subset, subset = nCount_RNA < 9e6 ) ##according to ???###
  
  out_file_data_subset=out_file_data[which(names(out_file_data) %in% colnames(RAW_subset))]
  RAW_subset=RAW_subset[,intersect(colnames(RAW_subset),names(out_file_data_subset))]
}  
  
tmp = RAW_subset@meta.data
data.dynamics = tmp[,-c(1:5,ncol(tmp))]
time_point
par(mfrow = c(3,3))
for(i in 1:9){
  timepoints = seq(0,time_point+5,5)
  plot(timepoints,data.dynamics[i,1:length(timepoints)],xlab = "time (min)",ylab = "Rela activity")
}

# plot rela activity dynamics
if(T){
  tmp = as.matrix(RAW_seurat@meta.data)
  TF.dyn = as.matrix(tmp[,-c(1:5,ncol(tmp))])
  TF.dyn = apply(TF.dyn,2, as.numeric)
  rownames(TF.dyn) = rownames(tmp)
  activity.m = apply(TF.dyn,2,function(x){mean(na.omit(x))})
  activity.sd = apply(TF.dyn,2,function(x){sd(na.omit(x))})
  thr = 0.2
  activity.qt = apply(TF.dyn,2,function(x){quantile(na.omit(x),c(thr,1-thr))})
  
  plot(seq(0,305,5),activity.m,type = "b",xlab = "time (min)",ylab = "Rela activity",ylim = c(0,1), main = "Rela activity dynamics")
  points(seq(0,305,5),activity.qt[1,],type = "b",col = "grey")
  points(seq(0,305,5),activity.qt[2,],type = "b",col = "grey")
  plot(seq(0,305,5),activity.sd/activity.m,type = "b",xlab = "time (min)",ylab = "Rela activity CV",ylim = c(0,1), main = "Rela activity CV")
}
# number of cells for each time point
apply(TF.dyn,2,function(x){return(length(na.omit(x)))})

# correlation between rela activity and TFA
data.dynamics = data.dynamics[order(rownames(data.dynamics)),]
ids = match(rownames(data.dynamics),rownames(df.all))
TFA.sub = df.all[ids,]

final_activity = data.dynamics[max_movie_time_point+1]
TFA.eval = cbind(TFA.sub[,1:2],final_activity)
colnames(TFA.eval)[3] = "TFA.gt"
cor(TFA.eval)

par(mfrow = c(1,1))
plot(as.vector(as.matrix(final_activity)),TFA.sub$TFA.in,xlab = "Rela activity",ylab = "TFA",col=  "red",ylim = c(0,10))
points(as.vector(as.matrix(final_activity)),TFA.sub$TFA.ex,col=  "blue")

# lm.fit = lm(as.vector(as.matrix(final_activity))~TFA.sub$TFA.in)
# summary(lm.fit)

# TFA.eval = cbind(TFA.sub[,1:2],apply(data.dynamics[1:max_movie_time_point+1],1,mean))
# colnames(TFA.eval)[3] = "TFA.gt"
# cor(TFA.eval)

# calculate Rela activity distribution in 4 timepoints
timepoint_Rela = c(1,max_movie_time_point_vector)

Rela_data = TF.dyn[,timepoint_Rela]
df.rela = data.frame(activity = as.vector(Rela_data),time = c(rep(0,nrow(TF.dyn)),rep(75,nrow(TF.dyn)),rep(150,nrow(TF.dyn)),
                                                              rep(300,nrow(TF.dyn))))
df.rela = na.omit(df.rela)
df.rela$time = factor(df.rela$time, levels = c(0,75,150,300))

final_activity_all = NULL
for(i in 1:nrow(TF.dyn)){
  tmp = na.omit(TF.dyn[i,])
  tmp = tmp[length(tmp)]
  final_activity_all = c(final_activity_all,tmp)
}
df.all[,6] = final_activity_all[match(rownames(df.all),rownames(TF.dyn))]
colnames(df.all)[6] = "Rela_activity"


# plot
library(ggplot2)
p1 = ggplot(df.all,aes(x = time,y = TFE)) + geom_boxplot(aes(fill = time)) 
p2 = ggplot(df.all,aes(x = time,y = TFA.in))+geom_boxplot(aes(fill = time)) 
p3 = ggplot(df.all,aes(x = time,y = TFA.ex))+geom_boxplot(aes(fill = time)) 
p4 = ggplot(df.rela,aes(x = time,y = activity))+geom_boxplot(aes(fill = time)) + ylim(0, 1)

library("cowplot")
pdf(paste(result_dir,"NFKB TFA dynamics boxplot",".pdf",sep = ""),height=9,width=9)
plot_grid(p4,p2,p3,p1,ncol = 1,align = "v")
dev.off()

# plot lines
df.rela.plot = data.frame(time = seq(0,300,5),mean = activity.m[-length(activity.m)],sd = activity.sd[-length(activity.m)])

p5 = ggplot(df.rela.plot,aes(x = time, y = mean)) + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),width = 0.2,
                                                                  position = position_dodge(0.1),colour = "purple")
p5 = p5 + geom_line(position = position_dodge(0.1),colour = "purple") + labs(y = "Nuclear Signal",x = "Time (min)")
pdf(paste(result_dir,"NFKB nuclear dynamics",".pdf",sep = ""),height=9,width=9)
print(p5)
dev.off()

# show only four time points
df.rela.plot.sub = df.rela.plot[which(is.element(df.rela.plot$time,c(0,75,150,300))),]
p1 = ggplot(df.rela.plot.sub,aes(x = time, y = mean)) + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),width = 0.2,
                                                                  position = position_dodge(0.1),colour = "purple")
p1 = p1 + geom_line(position = position_dodge(0.1),colour = "purple") + labs(y = "Nuclear Signal",x = NULL)

# calculate mean and sd
df.all.plot = NULL
for(i in 1:4){
  time_set = c(0,75,150,300)
  time = time_set[i]
  id = which(df.all$time==time)
  tmp = c(mean(df.all$TFA.in[id]),mean(df.all$TFA.ex[id]),mean(df.all$TFE[id],na.rm = T),mean(df.all$TFA.to[id],na.rm = T),
          sd(df.all$TFA.in[id]),sd(df.all$TFA.ex[id]),sd(df.all$TFE[id],na.rm = T),sd(df.all$TFA.to[id],na.rm = T))
  df.all.plot = rbind(df.all.plot,tmp)
}
colnames(df.all.plot) = c("intron.m","exon.m","exp.m","total.m","intron.sd","exon.sd","exp.sd","total.sd")
df.all.plot = as.data.frame(df.all.plot)
df.all.plot$time = time_set
# plot
p2 = ggplot(df.all.plot,aes(x = time, y = intron.m)) + geom_errorbar(aes(ymin = intron.m - intron.sd, ymax = intron.m + intron.sd),
                                                                width = 0.2, position = position_dodge(0.1),colour = "red")
p2 = p2 + geom_line(position = position_dodge(0.1),colour = "red")+labs(y = "Intron TFA",x = NULL)
p3 = ggplot(df.all.plot,aes(x = time, y = exon.m)) + geom_errorbar(aes(ymin = exon.m - exon.sd, ymax = exon.m + exon.sd),
                                                                     width = 0.2, position = position_dodge(0.1),colour = "blue")
p3 = p3 + geom_line(position = position_dodge(0.1),colour = "blue")+labs(y = "Exon TFA",x = NULL)
p4 = ggplot(df.all.plot,aes(x = time, y = exp.m)) + geom_errorbar(aes(ymin = exp.m - exp.sd, ymax = exp.m + exp.sd),
                                                                     width = 0.2, position = position_dodge(0.1),colour = "black")
p4 = p4 + geom_line(position = position_dodge(0.1),colour = "black")+labs(y = "TF Expression",x = "Time (min)")

library(cowplot)
pdf(paste(result_dir,"NFKB TFA dynamics",".pdf",sep = ""),height=9,width=9)
plot_grid(p1,p2,p3,p4,ncol = 1,align = "v")
dev.off()

colnames(df.rela.plot.sub) = c("time","nuclear.mean","nuclear.sd")
df.all.plot.ex = cbind(df.rela.plot.sub[,1:3],df.all.plot[,c(1,5,2,6,3,7,4,8)])
write.csv(df.all.plot.ex,file = paste(result_dir,"NFKB TFA dynamics",".csv",sep = ""),quote = F,row.names = F)

# compare intron, exon with intron+exon
if(F){
  # show intron TFA, exon TFA, total TFA in one cell
  df.all = log2(df.TFA+1)
  pdf(paste(result_dir,"intron+exon_CPM",".pdf",sep = ""),height=9,width=9)
  barplot(c(TFA_cells[1,1],TFA_cells[1,2],TFA_cells[1,1]+TFA_cells[1,2]),main = "Compare 3 types of TFA: NFkB",
          names.arg = c("intron_TFA","exon_TFA","total_TFA"), col = c("red","blue","green"),ylab = "CPM")
  dev.off()
  
  TFA_cells$TFA_totalCPM = TFA_cells$TFA_IntronCPM+TFA_cells$TFA_ExonCPM
  
  
  p5 = ggplot(df.all.plot,aes(x = time, y = total.m)) + geom_errorbar(aes(ymin = total.m - total.sd, ymax = total.m + total.sd),
                                                                      width = 0.2, position = position_dodge(0.1),colour = "green")
  p5 = p5 + geom_line(position = position_dodge(0.1),colour = "green")+labs(y = "Total TFA",x = NULL)
  pdf(paste(result_dir,"intron+exon_dynamics",".pdf",sep = ""),height=9,width=9)
  plot_grid(p1,p2,p3,p5,ncol = 1,align = "v")
  
  dev.off()
  
  
}

# show dynamics of individual genes
if(F){
  # load intronCPM and exonCPM data for rela targets
  TF_targets = colnames(mouse_dorothea_logic_matrix_ensembl)[which(mouse_dorothea_logic_matrix_ensembl["Rela",]==1)]
  if(F){ # ~5 min
    # record intron/exon CPM data
    out_file_dir="../raw data/tpmcalculator_out_file/"
    out_file_name_vector=list.files(out_file_dir)
    out_file_data=list()
    
    intronCPM_mat = exonCPM_mat = matrix(0,nrow = length(TF_targets),ncol = length(out_file_name_vector))
    rownames(intronCPM_mat) = rownames(exonCPM_mat) = TF_targets
    for(i in 1:length(out_file_name_vector)){
      out_file_i = out_file_name_vector[i]
      data_i=read.table(paste(out_file_dir,out_file_i,sep=""),fill = T,header=T)
      data_i=na.omit(data_i)
      rownames(data_i)=sub("\\..*$","",data_i$Gene_Id)
      data_i$ExonCPM=1e6*(data_i$ExonReads)/sum(data_i$Reads)
      data_i$IntronCPM=1e6*(data_i$IntronReads)/sum(data_i$Reads)
      
      over_genes = intersect(rownames(data_i),rownames(intronCPM_mat))
      exonCPM_mat[over_genes,i] = data_i[over_genes,"ExonCPM"]
      intronCPM_mat[over_genes,i] = data_i[over_genes,"IntronCPM"]
      
      
      print(i/length(out_file_name_vector))
    }
    # change names
    cell_sra_result=read.csv("../raw data/sra_result.csv")
    cell_SraRunInfo=read.csv("../raw data/SraRunInfo.csv")
    cell_metadata=merge(cell_sra_result,cell_SraRunInfo,all=F,by.x="Experiment.Accession",by.y="Experiment")
    rownames(cell_metadata)=cell_metadata$Run
    
    names_out_file_data_new=sub("_genes.out","",out_file_name_vector)
    names_out_file_data_new=cell_metadata[names_out_file_data_new,"Experiment.Title"]
    names_out_file_data_new=sub("\\).*$","",sub("^.*\\(","",names_out_file_data_new))
    names_out_file_data_new=paste("X",names_out_file_data_new,sep="")
    
    colnames(exonCPM_mat) = colnames(intronCPM_mat) = names_out_file_data_new
    # save(exonCPM_mat,intronCPM_mat,file = paste(result_dir,"intron_exon_CPM.RData",sep = ""))
    
  }
  load(paste(result_dir,"intron_exon_CPM.RData",sep = ""))
  # remove zero expression genes
  table(rowSums(exonCPM_mat)>0)
  table(rowSums(intronCPM_mat)>0)
  colnames(exonCPM_mat)
  # mean expression
  df.mean = data.frame(intron.m = apply(intronCPM_mat,1,mean),exon.m = apply(exonCPM_mat,1,mean))

  # metadata
  tp = c(0,75,150,300)
  time_data = RNA_cell_id$Time.point[na.omit(match(colnames(exonCPM_mat),rownames(RNA_cell_id)))]
  Rela_activity = final_activity_all[match(colnames(exonCPM_mat),rownames(TF.dyn))]
  
  dim(exonCPM_mat)
  
  # calculate correlation for each gene
  cor.genes = data.frame()
  for(gene_i in 1:nrow(exonCPM_mat)){
    data.t = data.frame()
    for(i in 1:4){
      time_set = c(0,75,150,300)
      time = time_set[i]
      id = which(time_data==time)
      tmp = c(mean(intronCPM_mat[gene_i,id]),mean(exonCPM_mat[gene_i,id]))
      data.t = rbind(data.t,tmp)
    }
    colnames(data.t) = c("intron.m","exon.m")
    # correlation with nuclear localization signal
    rela.t = df.rela.plot.sub$mean
    cor.in = cor(data.t$intron.m,rela.t)
    cor.ex = cor(data.t$exon.m,rela.t)
    
    cor.genes = rbind(cor.genes,c(cor.in,cor.ex))
  }
  rownames(cor.genes) = rownames(exonCPM_mat)
  colnames(cor.genes) = c("cor.in","cor.ex")
  
  # combine data
  df.comb = cbind(cor.genes,df.mean)
  # symbol names
  sym_names = geneIDs$SYMBOL[match(rownames(df.comb),geneIDs$GENEID)]
  df.comb$symbol = sym_names
  df.comb$symbol_toupper = toupper(sym_names)
  
  # evaluate the effect of mRNA half-life
  #import MCF7 half lives data
  MCF7_half_life_data=read.table("../../simulation data/raw code/p53 target half-lives/GSE49831_MCF7_halflives.txt",header=T)
  MCF7_half_life_data$MCF7_half_life_mean=apply(MCF7_half_life_data,1,mean)
  over_targets =intersect(toupper(sym_names),rownames(MCF7_half_life_data))
  mhl_sub = MCF7_half_life_data[over_targets,]$MCF7_half_life_mean
  names(mhl_sub) = over_targets
  
  # highest mRNA half-life and loweset mRNA half-life
  sort(mhl_sub)
  df.comb = cbind(df.comb,mhl = mhl_sub[df.comb$symbol_toupper])
  
  
  # expression distribution
  hist(log2(df.comb$intron.m+1))
  hist(log2(df.comb$exon.m+1))
  
  # select high expression genes
  sub.genes = which(df.comb$intron.m>1&df.comb$exon.m>1)
  df.sub = df.comb[sub.genes,]
  pdf(paste(result_dir,"indiv_correlation_nuclear",".pdf",sep = ""),height=5,width=5)
  
  plot(df.sub$cor.ex,df.sub$cor.in,xlab = "exonCPM correlation",ylab = "intronCPM correlation", main = "correlation with nuclear signal")
  abline(h=0,col = "red",lty = 2)
  abline(v=0,col = "blue",lty = 2)
  
  dev.off()
  # abline(a = 0,b = 1,col = "green",lty = 2)
  table((df.sub$cor.in)>(df.sub$cor.ex))
  table(abs(df.sub$cor.in)>abs(df.sub$cor.ex))
  table((df.sub$cor.in)>0 & (df.sub$cor.ex)>0)
  
  # bisection with mRNA half-life
  id.L = which(df.sub$mhl<median(df.sub$mhl,na.rm = T))
  id.H = which(df.sub$mhl>median(df.sub$mhl,na.rm = T))
  pdf(paste(result_dir,"indiv_boxplot",".pdf",sep = ""),height=5,width=5)
  
  boxplot(list(a = df.sub$cor.in[id.L], b = df.sub$cor.ex[id.L], c = df.sub$cor.in[id.H], d = df.sub$cor.ex[id.H]),
          names = c("L.in","L.ex","H.in","H.ex"),col = c("red","blue","red","blue"),main = "Effect of mRNA half-life",
          ylab = "correlation")
  
  dev.off()
  p1 = wilcox.test(df.sub$cor.in[id.L],df.sub$cor.ex[id.L],alternative = "greater")$p.value
  p2 = wilcox.test(df.sub$cor.in[id.H],df.sub$cor.ex[id.H],alternative = "greater")$p.value

  # pair-wise test
  id = which(df.sub$cor.in>0&df.sub$cor.ex>0&df.sub$mhl>0)
  df.sub2 = df.sub[id,]
  
  id.L = which(df.sub2$mhl<median(df.sub2$mhl,na.rm = T))
  id.H = which(df.sub2$mhl>median(df.sub2$mhl,na.rm = T))
  boxplot(list(a = df.sub2$cor.in[id.L]/df.sub2$cor.ex[id.L], b = df.sub2$cor.in[id.H]/df.sub2$cor.ex[id.H]),
          names = c("L","H"),main = "Effect of mRNA half-life",ylab = "correlation ratio: intron/exon")
  p1 = wilcox.test(df.sub2$cor.in[id.L]/df.sub2$cor.ex[id.L],
                   df.sub2$cor.in[id.H]/df.sub2$cor.ex[id.H],alternative = "less")$p.value
  p2 = t.test(df.sub2$cor.in[id.L]/df.sub2$cor.ex[id.L],
                   df.sub2$cor.in[id.H]/df.sub2$cor.ex[id.H],alternative = "less")$p.value
  
  
  
  # divide into 3 parts
  id.L = which(df.sub$mhl<quantile(df.sub$mhl,1/3,na.rm = T))
  id.M = which(df.sub$mhl<quantile(df.sub$mhl,2/3,na.rm = T)&df.sub$mhl>quantile(df.sub$mhl,1/3,na.rm = T))
  id.H = which(df.sub$mhl>quantile(df.sub$mhl,2/3,na.rm = T))
  boxplot(list(a = df.sub$cor.in[id.L], b = df.sub$cor.ex[id.L], c = df.sub$cor.in[id.M], d = df.sub$cor.ex[id.M],
               e = df.sub$cor.in[id.H], f = df.sub$cor.ex[id.H]),
          names = c("L.in","L.ex","M.in","M.ex","H.in","H.ex"),
          col = c("red","blue","red","blue","red","blue"),main = "Effect of mRNA half-life",
          ylab = "correlation")
  
  p1 = wilcox.test(df.sub$cor.in[id.L],df.sub$cor.ex[id.L],alternative = "greater")$p.value
  p2 = wilcox.test(df.sub$cor.in[id.H],df.sub$cor.ex[id.H],alternative = "greater")$p.value
  
  # calculate correlation between mRNA half-life and intron/exon performance
  plot(df.sub$mhl,df.sub$cor.in/df.sub$cor.ex)
  # select only positive correlation
  id = which(df.sub$cor.in>0&df.sub$cor.ex>0&df.sub$mhl>0)
  df.sub2 = df.sub[id,]
  pdf(paste(result_dir,"indiv_mhl_cor",".pdf",sep = ""),height=5,width=5)
  
  plot(df.sub2$mhl,df.sub2$cor.in/df.sub2$cor.ex,xlab = "mRNA half-life (min)", ylab = "cor.in/cor.ex",log ="x",
       main = "Effect of mRNA half-life")
  
  dev.off()
  cor(df.sub2$mhl,df.sub2$cor.in/df.sub2$cor.ex)
  
  cor.test(df.sub2$mhl,df.sub2$cor.in/df.sub2$cor.ex)
  gene1 = rownames(df.sub)[which(df.sub$cor.in==max(df.sub$cor.in))]
  gene2 = rownames(df.sub)[which(df.sub$cor.ex==max(df.sub$cor.ex))]
  gene3 = rownames(df.sub)[which(df.sub$cor.in==min(df.sub$cor.in))]
  
  # low mRNA half-life gene and high mRNA half-life gene
  gene4 = rownames(df.sub[order(df.sub$mhl,decreasing = F),])[1]
  gene5 = rownames(df.sub[order(df.sub$mhl,decreasing = T),])[1]
  
  
  # visualize one gene
  gene_i = gene5
  if(T){
    # calculate mean and sd
    df.all.plot = NULL
    for(i in 1:4){
      time_set = c(0,75,150,300)
      time = time_set[i]
      id = which(time_data==time)
      tmp = c(mean(intronCPM_mat[gene_i,id]),mean(exonCPM_mat[gene_i,id]),
              sd(intronCPM_mat[gene_i,id]),sd(exonCPM_mat[gene_i,id]))
      df.all.plot = rbind(df.all.plot,tmp)
    }
    colnames(df.all.plot) = c("intron.m","exon.m","intron.sd","exon.sd")
    df.all.plot = as.data.frame(df.all.plot)
    df.all.plot$time = time_set
    
    p1 = ggplot(df.rela.plot.sub,aes(x = time, y = mean)) + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),width = 0.2,
                                                                          position = position_dodge(0.1),colour = "purple")
    p1 = p1 + geom_line(position = position_dodge(0.1),colour = "purple") + 
      labs(y = "Nuclear Signal",x = NULL,title = df.sub[gene_i,]$symbol)
    
    
    p2 = ggplot(df.all.plot,aes(x = time, y = intron.m)) + geom_errorbar(aes(ymin = intron.m - intron.sd, ymax = intron.m + intron.sd),
                                                                         width = 0.2, position = position_dodge(0.1),colour = "red")
    p2 = p2 + geom_line(position = position_dodge(0.1),colour = "red")+labs(y = "IntronCPM",x = NULL)
    p3 = ggplot(df.all.plot,aes(x = time, y = exon.m)) + geom_errorbar(aes(ymin = exon.m - exon.sd, ymax = exon.m + exon.sd),
                                                                       width = 0.2, position = position_dodge(0.1),colour = "blue")
    p3 = p3 + geom_line(position = position_dodge(0.1),colour = "blue")+labs(y = "ExonCPM",x = NULL)
    plot_grid(p1,p2,p3,ncol = 1,align = "v")
    
  }
  # visulize individual cells
  time = 75
  id = which(time_data==time)
  mat = intronCPM_mat[,id]
  # select gens with highest intron expression
  mat.m = sort(apply(mat,1,mean),decreasing =T) 
  sel_genes = names(mat.m)[1:9]
  
  # Rela activity data
  if(T){
    time_point_vector=c(75,150,300)
    max_movie_time_point_vector=time_point_vector/5+1
    i_time_point = 1
    time_point=time_point_vector[i_time_point]
    max_movie_time_point=max_movie_time_point_vector[i_time_point]
    RAW_subset=subset(RAW_seurat, subset = Time.point == time_point)
    RAW_subset=subset(RAW_subset, subset = nCount_RNA < 9e6 ) ##according to ???###
    
    out_file_data_subset=out_file_data[which(names(out_file_data) %in% colnames(RAW_subset))]
    RAW_subset=RAW_subset[,intersect(colnames(RAW_subset),names(out_file_data_subset))]
    
    tmp = RAW_subset@meta.data
    data.dynamics = tmp[,-c(1:5,ncol(tmp))]
    time_point
    par(mfrow = c(3,3))
    for(i in 1:9){
      timepoints = seq(0,time_point+5,5)
      plot(timepoints,data.dynamics[i,1:length(timepoints)],xlab = "time (min)",ylab = "Rela activity")
    }
    
    # number of cells for each time point
    apply(TF.dyn,2,function(x){return(length(na.omit(x)))})
    
    # correlation between rela activity and TFA
    data.dynamics = data.dynamics[order(rownames(data.dynamics)),]
    ids = match(rownames(data.dynamics),rownames(df.all))
    TFA.sub = df.all[ids,]
    
    final_activity = data.dynamics[max_movie_time_point+1]
    TFA.eval = cbind(TFA.sub[,1:2],final_activity)
    colnames(TFA.eval)[3] = "TFA.gt"
    cor(TFA.eval)
    
    par(mfrow = c(1,1))
    plot(as.vector(as.matrix(final_activity)),TFA.sub$TFA.in,xlab = "Rela activity",ylab = "TFA",col=  "red",ylim = c(0,10))
    points(as.vector(as.matrix(final_activity)),TFA.sub$TFA.ex,col=  "blue")
    cor(as.vector(as.matrix(final_activity)),TFA.sub$TFA.in)
    cor(as.vector(as.matrix(final_activity)),TFA.sub$TFA.ex)
    
  }
  
  # calculate TFA for high half-life and low half-life genes
  if(T){
    time = 75
    id = which(time_data==time)
    mat = intronCPM_mat[,id]
    
    id.tmp = which(!is.na(df.comb$mhl))
    df.tmp = df.comb[id.tmp,]
    genes.H = rownames(df.tmp)[which(df.tmp$mhl>median(df.tmp$mhl))]
    genes.L = rownames(df.tmp)[which(df.tmp$mhl<median(df.tmp$mhl))]
    TFA_in.H = log2(apply(intronCPM_mat[genes.H,id],2,mean)+1)
    TFA_ex.H = log2(apply(exonCPM_mat[genes.H,id],2,mean)+1)
    TFA_in.L = log2(apply(intronCPM_mat[genes.L,id],2,mean)+1)
    TFA_ex.L = log2(apply(exonCPM_mat[genes.L,id],2,mean)+1)
    
    pdf(paste(result_dir,"TFA dynamics and correlation",".pdf",sep = ""),height=5,width=5)
    
    plot(as.vector(as.matrix(final_activity)),TFA_in.H,xlab = "Rela activity",ylab = "TFA",col=  "red",ylim = c(0,10),
         main = "High mRNA half-life")
    points(as.vector(as.matrix(final_activity)),TFA_ex.H,col=  "blue")
    
    plot(as.vector(as.matrix(final_activity)),TFA_in.L,xlab = "Rela activity",ylab = "TFA",col=  "red",ylim = c(0,10),
         main = "Low mRNA half-life")
    points(as.vector(as.matrix(final_activity)),TFA_ex.L,col=  "blue")
    
    cors = c(cor(as.vector(as.matrix(final_activity)),TFA_in.H),cor(as.vector(as.matrix(final_activity)),TFA_ex.H),
             cor(as.vector(as.matrix(final_activity)),TFA_in.L),cor(as.vector(as.matrix(final_activity)),TFA_ex.L))
    names(cors) = c("H","H","L","L")
    cols = c("red","blue","red","blue")
    barplot(cors,col = cols,main = "correlation with Rela signal")
    
    dev.off()
  }
  
  # show scatter plot of gene dynamics with TF dynamics
  pdf(paste(result_dir,"sc_dynamics",".pdf",sep = ""),height=5,width=5)
  
  par(mfrow = c(1,1))
  df.cor = NULL
  for(i in 1:length(sel_genes)){
    gene = sel_genes[i]
    ids = match(rownames(data.dynamics),colnames(intronCPM_mat))
    vec.in = intronCPM_mat[gene,ids]
    vec.ex = exonCPM_mat[gene,ids]
    
    # symbol names
    sym_name = geneIDs$SYMBOL[match(gene,geneIDs$GENEID)]
    
    plot(as.vector(as.matrix(final_activity)),vec.in,xlab = "Rela activity",ylab = "CPM",col=  "red",main = sym_name,
         ylim= range(c(vec.in,vec.ex)))
    points(as.vector(as.matrix(final_activity)),vec.ex,col=  "blue")
    cor1 = cor.test(as.vector(as.matrix(final_activity)),vec.in)
    cor2 = cor.test(as.vector(as.matrix(final_activity)),vec.ex)
    
    df.cor = rbind(df.cor,c(cor1$estimate,cor1$p.value,cor2$estimate,cor2$p.value))
  }
  colnames(df.cor) = c("cor.in","pvalue.in","cor.ex","pvalue.ex")
  rownames(df.cor) = geneIDs$SYMBOL[match(sel_genes,geneIDs$GENEID)]
  df.cor = data.frame(df.cor)
  dev.off()
}
