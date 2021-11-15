rm(list = ls())
start_time=Sys.time()
library("pheatmap")
library("reshape2")
library("ggplot2")
library("readxl")
## install multiplot
# library(devtools) 
# install_github( "kassambara/easyGgplot2")
library("easyGgplot2")
source("../../materials//scTFA_calculation.R")

themo_demo=theme(
  #plot.margin = margin(2, 2, 2, 2, "cm"),
  text = element_text(size=30),
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_blank(),
  axis.text=element_text(color='black'),
  plot.title = element_text(hjust = 0.5),
  title=element_text(size = 30),
  axis.title.x = element_text(size=24),
  axis.title.y = element_text(size=24),
  axis.text.x = element_text(size=22,angle=45,hjust = 1),
  axis.text.y = element_text(size=22))

data_dir="../raw data/"
result_dir="../results/"
date_analysis="analyzed20210129_1(scTFA)"
dir.create(paste(result_dir))

#import GRN
human_dorothea_logic_matrix=read.csv("../../materials/human_GRN_classAB_activation_only/dorothea_human_AB.csv",
                                     row.names = 1,check.names=FALSE)

#import metadata
sra_result=read.csv("../raw data/sra_result.csv")
SraRunInfo=read.csv("../raw data/SraRunInfo.csv")
  
#import tpmcalculator data
out_file_name_vector=list.files(paste(data_dir,"tpmcalculator_result_analyzed20200618_1",sep = ""))
out_file_name_vector=grep("genes.out$",out_file_name_vector,value = T)
index=0
for(out_file_i in out_file_name_vector){
  data_i=read.table(paste(data_dir,"tpmcalculator_result_analyzed20200618_1/",out_file_i,sep=""),fill = T,header=T)
  data_i=na.omit(data_i)
  data_i=data_i[!duplicated(data_i$Gene_Id),]
  rownames(data_i)=data_i$Gene_Id

  data_i$ExonCPM=1e6*(data_i$ExonReads)/sum(data_i$Reads)
  data_i$IntronCPM=1e6*(data_i$IntronReads)/sum(data_i$Reads)
  
  data_i_chosen=data_i[,c("Gene_Id","IntronLength","ExonCPM","IntronCPM")]
  
  #add cell name
  SRR_name=sub("_Aligned.sortedByCoord.out_genes.out","",out_file_i)
  SRR_name=sub("\\..*$","",SRR_name)
  SRX_name=SraRunInfo[which(SraRunInfo$Run==SRR_name),"Experiment"]
  sample_name=sra_result[which(sra_result$Experiment.Accession==SRX_name),"Experiment.Title"]
  sample_name=sub("; Homo sapiens; RNA-Seq","",sample_name)
  sample_name=sub("^.*RNA-Seq, ","",sample_name)
  colnames(data_i_chosen)[-1]=paste(colnames(data_i_chosen)[-1],sample_name,sep="_")
  
  if(index==0){
    data_all_combined=data_i_chosen
  }else{
    data_all_combined=merge(data_all_combined,data_i_chosen,by.x="Gene_Id",by.y="Gene_Id",all=T)
  }
  index=index+1
  
  print(index/length(out_file_name_vector))
}
rownames(data_all_combined)=data_all_combined$Gene_Id
write.csv(data_all_combined,paste(result_dir,"data_all_combined",date_analysis,".csv"))

# chosen genes in GRN
data_all_combined_GRN=data_all_combined[data_all_combined$Gene_Id %in% colnames(human_dorothea_logic_matrix),]
# preserve genes with non-NA value in all samples, then change NA to 0
data_all_combined_GRN_subset=data_all_combined_GRN[rowSums(!is.na(data_all_combined_GRN))>=
                                                     ncol(data_all_combined_GRN),]
data_all_combined_GRN_subset[is.na(data_all_combined_GRN_subset)]=0

IntronLength_vector=apply(data_all_combined_GRN_subset[,grep("IntronLength",colnames(data_all_combined_GRN_subset))],
                          1,function(x) mean(x[!is.na(x)]))

data_IntronCPM=data_all_combined_GRN_subset[IntronLength_vector>0,grep("IntronCPM",colnames(data_all_combined_GRN_subset),value = T)]
data_IntronCPM=data_IntronCPM[,-(grep("t=24 h",colnames(data_IntronCPM)))] #remove 24h data


TFA_IntronCPM=scTFA_calculation(as.matrix(data_IntronCPM),human_dorothea_logic_matrix,zscore=F,
                                                   gene_weight_method=NULL)
data_ExonCPM=data_all_combined_GRN_subset[,grep("ExonCPM",colnames(data_all_combined_GRN_subset),value = T)]
data_ExonCPM=data_ExonCPM[,-(grep("t=24 h",colnames(data_ExonCPM)))] #remove 24h data
TFA_ExonCPM=scTFA_calculation(as.matrix(data_ExonCPM),human_dorothea_logic_matrix,zscore=F,
                                                   gene_weight_method=NULL)

write.csv(data_IntronCPM,paste(result_dir,"data_IntronCPM",date_analysis,".csv"))
write.csv(data_ExonCPM,paste(result_dir,"data_ExonCPM",date_analysis,".csv"))

write.csv(TFA_IntronCPM,paste(result_dir,"TFA_IntronCPM",date_analysis,".csv"))
write.csv(TFA_ExonCPM,paste(result_dir,"TFA_ExonCPM",date_analysis,".csv"))

TFA_all_combined=cbind(TFA_IntronCPM,TFA_ExonCPM[rownames(TFA_IntronCPM),])
write.csv(TFA_all_combined,paste(result_dir,"TFA_all_combined",date_analysis,".csv"))

rep=c("rep1","rep2")
TF_chosen="TP53"
data_chosen=as.data.frame(TFA_all_combined["TP53",])
colnames(data_chosen)=TF_chosen
data_chosen$index=sub("_t=.*$","",rownames(data_chosen))
data_chosen$normalization=sub("Intron","",sub("Exon","",data_chosen$index))
data_chosen$rep=sub("_1","",sub("^.*, ","",rownames(data_chosen)))
data_chosen$condition=sub("\\, rep.*$","",sub("^.*t=","",rownames(data_chosen)))
data_chosen$condition=factor(data_chosen$condition,levels=c("0 h",paste(c(1:12)," h, IR 10Gy",sep="")))
write.csv(data_chosen,paste(result_dir,"data_chosen",TF_chosen,date_analysis,".csv"))

# target number of P53
rowSums(human_dorothea_logic_matrix)["TP53"]
gt_targets = colnames(human_dorothea_logic_matrix)[which(human_dorothea_logic_matrix["TP53",]==1)]
targets_in_exon_data = intersect(gt_targets,rownames(data_ExonCPM))
targets_in_intron_data = intersect(gt_targets,rownames(data_IntronCPM))

method_vector=c("ExonCPM","IntronCPM")
data_chosen_matrix=matrix(,length(levels(data_chosen$condition)),length(method_vector)*length(rep),
                          dimnames = list(levels(data_chosen$condition),
                                          paste(rep(method_vector,each=length(rep)),rep,sep="_")))

data_chosen_matrix_average=matrix(,length(levels(data_chosen$condition)),length(method_vector),
                                  dimnames = list(levels(data_chosen$condition),paste(method_vector,"average",sep="_")))

for(row_i in rownames(data_chosen_matrix)){
  for(method_j in method_vector){
    for(rep_k in rep){
      data_chosen_matrix[row_i,paste(method_j,rep_k,sep="_")]=
        data_chosen[which((data_chosen$condition==row_i)*(data_chosen$index==method_j)*(data_chosen$rep==rep_k)==1),TF_chosen]
    }
    data_chosen_matrix_average[row_i,paste(method_j,"average",sep="_")]=
      mean(data_chosen[which((data_chosen$condition==row_i)*(data_chosen$index==method_j)==1),TF_chosen])
  }
}

data_chosen_matrix=cbind(data_chosen_matrix,data_chosen_matrix_average[rownames(data_chosen_matrix),])


# import P53 ChIP-seq data
P53_chip_data=as.data.frame(read_xlsx("../raw data/41594_2017_BFnsmb3452_MOESM4_ESM.xlsx",sheet = 1))
# subset data to preserve peaks within 2kb of TSS
P53_chip_data_subset_2kb=P53_chip_data[abs(P53_chip_data$DisttoTSS)<=2000,]
P53_target=colnames(human_dorothea_logic_matrix)[human_dorothea_logic_matrix["TP53",]==1]
P53_target_in_chip=intersect(P53_target,P53_chip_data_subset_2kb$ClosestGene)
# chip binding closest genes
ChIP_genes = unique(P53_chip_data_subset_2kb$ClosestGene)

P53_chip_data_subset_2kb_target_gene=P53_chip_data_subset_2kb[P53_chip_data_subset_2kb$ClosestGene %in% P53_target_in_chip,]

col_chosen=c("NormReads_p53_t0","NormReads_p53_t1","NormReads_p53_t2_5",
             "NormReads_p53_t4","NormReads_p53_t5","NormReads_p53_t7_5")
P53_chip_data_summary=data.frame(chip_condition=col_chosen,
                                 chip_data_only_target_within_2kb=apply(P53_chip_data_subset_2kb_target_gene[,col_chosen],2,sum),
                                 chip_data_all=apply(P53_chip_data[,col_chosen],2,sum))
P53_chip_data_summary$chip_data_only_target_within_2kb_FC=
  P53_chip_data_summary$chip_data_only_target_within_2kb/
  P53_chip_data_summary$chip_data_only_target_within_2kb[1]
P53_chip_data_summary$chip_data_all_FC=
  P53_chip_data_summary$chip_data_all/
  P53_chip_data_summary$chip_data_all[1]

write.csv(P53_chip_data_summary,paste(result_dir,"P53_chip_data_summary",date_analysis,".csv"))

line_plot_data=data_chosen_matrix[,c("IntronCPM_average","ExonCPM_average")]
line_plot_data=as.data.frame(line_plot_data)
line_plot_data$condition=factor(rownames(line_plot_data),levels=rownames(line_plot_data))
line_plot_data_FC=line_plot_data
line_plot_data_FC$IntronCPM_average=line_plot_data_FC$IntronCPM_average/line_plot_data_FC$IntronCPM_average[1]
line_plot_data_FC$ExonCPM_average=line_plot_data_FC$ExonCPM_average/line_plot_data_FC$ExonCPM_average[1]

line_plot_data_FC_data_frame=melt(line_plot_data_FC,id.vars = "condition",variable.name = 'index')
line_plot_data_FC_data_frame$time=as.numeric(sub(" .*$","",as.vector(line_plot_data_FC_data_frame$condition)))

P53_chip_data_summary_plot=melt(P53_chip_data_summary[,c("chip_condition","chip_data_only_target_within_2kb_FC","chip_data_all_FC")],
                                id.vars = "chip_condition",variable.name = 'index')
P53_chip_data_summary_plot$time=as.numeric(sub("_",".",sub("NormReads_p53_t","",P53_chip_data_summary_plot$chip_condition)))

line_plot_data_FC_data_frame=rbind(line_plot_data_FC_data_frame[,c("time","index","value")],
                                P53_chip_data_summary_plot[,c("time","index","value")])

p1<-ggplot(data=line_plot_data_FC_data_frame,mapping=aes(x=time,y=value,fill=index,color=index,group=index))+
  geom_line(size=2)+
  geom_point(size = 5)+
  labs(x = "condition", y = "TFA_average ChIP_FC",
       title = TF_chosen)+
  theme_bw()+themo_demo

# linear Interpolation
time_point_final=seq(0,7.5,by=0.5)
index=0
for(index_i in names(table(line_plot_data_FC_data_frame$index))){
  data_index_i=line_plot_data_FC_data_frame[line_plot_data_FC_data_frame$index==index_i,]
  data_index_i=data_index_i[order(data_index_i$time),]
  data_index_i_interpolated=data_index_i[1,]
  for(time_point_i in time_point_final){
    if(time_point_i %in% data_index_i$time){
      data_index_i_interpolated=rbind(data_index_i_interpolated,data_index_i[data_index_i$time==time_point_i,])
    }else{
      time_point_before=data_index_i[max(which(data_index_i$time<time_point_i)),"time"]
      value_before=data_index_i[max(which(data_index_i$time<time_point_i)),"value"]
      
      time_point_after=data_index_i[min(which(data_index_i$time>time_point_i)),"time"]
      value_after=data_index_i[min(which(data_index_i$time>time_point_i)),"value"]
      
      value_now=value_before+((value_after-value_before)/(time_point_after-time_point_before))*(time_point_i-time_point_before)
      data_index_i_interpolated=rbind(data_index_i_interpolated,c(time_point_i,index_i,value_now))
    }
  }
  data_index_i_interpolated=data_index_i_interpolated[-1,]
  if(index==0){
    line_plot_data_FC_data_frame_interpolated=data_index_i_interpolated
  }else{
    line_plot_data_FC_data_frame_interpolated=rbind(line_plot_data_FC_data_frame_interpolated,
                                                 data_index_i_interpolated)
  }
  index=index+1
}
line_plot_data_FC_data_frame_interpolated$value=as.numeric(line_plot_data_FC_data_frame_interpolated$value)

p2<-ggplot(data=line_plot_data_FC_data_frame_interpolated,mapping=aes(x=time,y=value,fill=index,color=index,group=index))+
  geom_line(size=2)+
  geom_point(size = 5)+
  labs(x = "condition", y = "TFA_average ChIP_FC",
       title = TF_chosen)+
  theme_bw()+themo_demo
pdf(file =paste(result_dir,"TFA and ChIP_FC of",TF_chosen,"interpolated",date_analysis,".pdf"),width = 30,height = 10)
ggplot2.multiplot(p1,p2,cols=2)
dev.off()

line_plot_data_FC_data_frame_interpolated_dcast=dcast(line_plot_data_FC_data_frame_interpolated,index~time)
rownames(line_plot_data_FC_data_frame_interpolated_dcast)=line_plot_data_FC_data_frame_interpolated_dcast$index
line_plot_data_FC_data_frame_interpolated_dcast=t(line_plot_data_FC_data_frame_interpolated_dcast[,-1])

pdf(paste(result_dir,"heatmap of correlation between TFA and chip data of",TF_chosen,date_analysis,".pdf"),height=6,width=6)
pheatmap(cor(line_plot_data_FC_data_frame_interpolated_dcast),scale="none",
         cellwidth = 45, cellheight = 45,
         cluster_rows = F,cluster_cols = F,display_numbers = T,
         fontsize_number = 15,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(101),breaks=seq(0,1,length.out = 101),
         main=paste("Pearson correlation"), angle_col=90)
graphics.off()

# end_time=Sys.time()
# capture.output(end_time-start_time,file=paste("runtime",date_analysis,".txt"))
# capture.output(sessionInfo(),file=paste("sessionInfo",date_analysis,".txt"))

# compare TFA with TF expression

# expression of p53
TF_exp_rep1 = data_ExonCPM["TP53",grep("rep1",colnames(data_ExonCPM))]
TF_exp_rep2 = data_ExonCPM["TP53",grep("rep2",colnames(data_ExonCPM))]
TF_exp_rep1_FC = TF_exp_rep1/as.numeric(TF_exp_rep1[1])
TF_exp_rep2_FC = TF_exp_rep2/as.numeric(TF_exp_rep2[1])
plot(0:12,TF_exp_rep1_FC,col = "blue",type = "b",ylim = c(0.8,4))
points(0:12,TF_exp_rep2_FC,col = "green",type = "b")


# plot IntronTFA and ExonTFA in two replicates
TFA_in_rep1 = TFA_IntronCPM["TP53",grep("rep1",colnames(TFA_IntronCPM))]
TFA_in_rep2 = TFA_IntronCPM["TP53",grep("rep2",colnames(TFA_IntronCPM))]
TFA_in_rep1_FC = TFA_in_rep1/as.numeric(TFA_in_rep1[1])
TFA_in_rep2_FC = TFA_in_rep2/as.numeric(TFA_in_rep2[1])

TFA_ex_rep1 = TFA_ExonCPM["TP53",grep("rep1",colnames(TFA_ExonCPM))]
TFA_ex_rep2 = TFA_ExonCPM["TP53",grep("rep2",colnames(TFA_ExonCPM))]
TFA_ex_rep1_FC = TFA_ex_rep1/as.numeric(TFA_ex_rep1[1])
TFA_ex_rep2_FC = TFA_ex_rep2/as.numeric(TFA_ex_rep2[1])

t_chip = c(0,1,2.5,4,5,7.5)
chip_signals = P53_chip_data_summary$chip_data_all_FC

par(mfrow = c(1,2))
plot(0:12,TFA_in_rep1_FC,col = "red",type = "b",log = "y",ylim = c(0.8,5),main = "rep1 dynamics",xlab = "time (h)")
points(0:12,TFA_ex_rep1_FC,col = "blue",type = "b")
points(t_chip,chip_signals,col = "purple",type = "b")
points(0:12,TF_exp_rep1_FC,col = "black",type = "b",lty = 2)

plot(0:12,TFA_in_rep2_FC,col = "red",type = "b",log = "y",ylim = c(0.8,5),main = "rep2 dynamics",xlab = "time (h)")
points(0:12,TFA_ex_rep2_FC,col = "blue",type = "b")
points(t_chip,chip_signals,col = "purple",type = "b")
points(0:12,TF_exp_rep2_FC,col = "black",type = "b",lty = 2)

# compare expression level of p53
par(mfrow = c(1,1))
plot(0:12,TF_exp_rep1,col = "black",type = "b",main = "p53 expression dynamics (ExonCPM)",xlab = "time (h)",ylim = c(0,80))
points(0:12,TF_exp_rep2,col = "black",type = "b",pch = 2)


# Linear Interpolation
LinearInterpolation = function(t,x,t_final){
  x_final = NULL
  t = sort(t)
  for(i in 1:length(t_final)){
    if(is.element(t_final[i],t)){
      x_final = c(x_final,as.numeric(x[which(t==t_final[i])]))
    }else{
      id.max = which(t>t_final[i])[1]
      id.min = id.max - 1
      k = (x[id.max]-x[id.min])/(t[id.max]-t[id.min])
      x_final = c(x_final,as.numeric(k*(t_final[i]-t[id.max])+x[id.max]))
    }
  }
  return(x_final)
}
t_final = seq(0,7.5,0.5)
chip_LI = LinearInterpolation(t_chip,chip_signals,t_final)
TFA_in_rep1_LI = LinearInterpolation(0:12,TFA_in_rep1,t_final)
TFA_ex_rep1_LI = LinearInterpolation(0:12,TFA_ex_rep1,t_final)
TFE_rep1_LI = LinearInterpolation(0:12,TF_exp_rep1,t_final)
TFA_in_rep2_LI = LinearInterpolation(0:12,TFA_in_rep2,t_final)
TFA_ex_rep2_LI = LinearInterpolation(0:12,TFA_ex_rep2,t_final)
TFE_rep2_LI = LinearInterpolation(0:12,TF_exp_rep2,t_final)

TFA_in_rep1_LI_FC= TFA_in_rep1_LI/TFA_in_rep1_LI[1]
TFA_ex_rep1_LI_FC = TFA_ex_rep1_LI/TFA_ex_rep1_LI[1]
TFE_rep1_LI_FC = TFE_rep1_LI/TFE_rep1_LI[1]
TFA_in_rep2_LI_FC = TFA_in_rep2_LI/TFA_in_rep2_LI[1]
TFA_ex_rep2_LI_FC = TFA_ex_rep2_LI/TFA_ex_rep2_LI[1]
TFE_rep2_LI_FC = TFE_rep2_LI/TFE_rep2_LI[1]


par(mfrow = c(1,2))
plot(t_final,TFA_in_rep1_LI_FC,col = "red",type = "b",log = "y",ylim = c(0.8,5),main = "rep1 dynamics (interpolation)",xlab = "time (h)")
points(t_final,TFA_ex_rep1_LI_FC,col = "blue",type = "b")
points(t_final,chip_LI,col = "purple",type = "b")
points(t_final,TFE_rep1_LI_FC,col = "black",type = "b",lty = 2)

plot(t_final,TFA_in_rep2_LI_FC,col = "red",type = "b",log = "y",ylim = c(0.8,5),main = "rep2 dynamics (interpolation)",xlab = "time (h)")
points(t_final,TFA_ex_rep2_LI_FC,col = "blue",type = "b")
points(t_final,chip_LI,col = "purple",type = "b")
points(t_final,TFE_rep2_LI_FC,col = "black",type = "b",lty = 2)

df.rep1 = data.frame(TFA.in = TFA_in_rep1_LI_FC,TFA.ex = TFA_ex_rep1_LI_FC, TFE = TFE_rep1_LI_FC,chip = chip_LI)
cor(df.rep1)
df.rep2 = data.frame(TFA.in = TFA_in_rep2_LI_FC,TFA.ex = TFA_ex_rep2_LI_FC, TFE = TFE_rep2_LI_FC,chip = chip_LI)
cor(df.rep2)

# plot mean results
TFA_in_mean = (TFA_in_rep1_LI + TFA_in_rep2_LI)/2
TFA_in_mean = TFA_in_mean/TFA_in_mean[1]
TFA_ex_mean = (TFA_ex_rep1_LI + TFA_ex_rep2_LI)/2
TFA_ex_mean = TFA_ex_mean/TFA_ex_mean[1]
TFE_mean = (TFE_rep1_LI + TFE_rep2_LI)/2
TFE_mean = TFE_mean/TFE_mean[1]
par(mfrow = c(1,1))
pdf(paste(result_dir,"P53 TFA dynamics",".pdf",sep = ""),height=9,width=9)
plot(t_final,TFA_in_mean,col = "red",type = "b",log = "y",ylim = c(0.8,5),main = "P53 Dynamics"
     ,xlab = "Time (h)",ylab = "Fold Change")
points(t_final,TFA_ex_mean,col = "blue",type = "b")
points(t_final,chip_LI,col = "purple",type = "b")
points(t_final,TFE_mean,col = "black",type = "b")
dev.off()
df.mean = data.frame(TFA.in = TFA_in_mean,TFA.ex = TFA_ex_mean, TFE = TFE_mean,chip = chip_LI)
cor(df.mean)

# plot correlation results
df.cor = rbind(cor(df.rep1)["chip",1:3],cor(df.rep2)["chip",1:3],cor(df.mean)["chip",1:3])
rownames(df.cor) = c("rep1","rep2","mean")
df.plot = data.frame(cor = c(df.cor[1,],df.cor[2,],df.cor[3,]), sample = c(rep("rep1",3),rep("rep2",3),rep("mean",3)),
                     data.type = rep(c("TFA.intron","TFA.exon","TFexpression"),3))
df.plot$sample = factor(df.plot$sample,levels = c("rep1","rep2","mean"))
df.plot$data.type = factor(df.plot$data.type,levels = c("TFA.intron","TFA.exon","TFexpression"))

pdf(paste(result_dir,"Correlation with ChIP",".pdf",sep = ""),height=9,width=9)
library(ggplot2)
p = ggplot() + geom_bar(data = df.plot, aes(x = sample,y = cor,  fill = data.type),stat = "identity",
                        position = "dodge")+ scale_fill_manual(values=c("red","blue","black"))

p = p + ylim(-1,1) + labs(x = NULL,y = "correlation",title = "Correlation with ChIP data")
p + theme(plot.title = element_text(size = 20))
dev.off()

write.csv(df.plot,file = paste(result_dir,"Correlation with ChIP",".csv",sep = ""),quote = F,row.names = F)

# show total TFA results
TFA_to_mean = (TFA_ex_rep1_LI + TFA_ex_rep2_LI+TFA_in_rep1_LI + TFA_in_rep2_LI)/2
TFA_to_mean = TFA_to_mean/TFA_to_mean[1]
pdf(paste(result_dir,"P53 TFA dynamics",".pdf",sep = ""),height=9,width=9)
plot(t_final,TFA_in_mean,col = "red",type = "b",log = "y",ylim = c(0.8,5),main = "P53 Dynamics"
     ,xlab = "Time (h)",ylab = "Fold Change")
points(t_final,TFA_ex_mean,col = "blue",type = "b")
points(t_final,chip_LI,col = "purple",type = "b")
points(t_final,TFA_to_mean,col = "green",type = "b")

# calculate TFA with different set of genes
if(T){
  GRN=read.csv("../../materials/human_GRN_classAB_activation_only/dorothea_human_AB.csv",
               row.names = 1,check.names=FALSE)
  
  TF_targets = colnames(GRN)[which(GRN['TP53',]==1)]
  # evaluate the effect of mRNA half-life
  #import MCF7 half lives data
  MCF7_half_life_data=read.table("../../simulation data/raw code/p53 target half-lives/GSE49831_MCF7_halflives.txt",header=T)
  MCF7_half_life_data$MCF7_half_life_mean=apply(MCF7_half_life_data,1,mean)
  over_targets =intersect(TF_targets,rownames(MCF7_half_life_data))
  mhl_sub = MCF7_half_life_data[over_targets,]$MCF7_half_life_mean
  names(mhl_sub) = over_targets
  mhl_sub = sort(mhl_sub)
  plot(1:length(mhl_sub),mhl_sub,xlab = "id",ylab = "mRNA half-life (min)",log = 'y')
  gene.L = names(mhl_sub)[which(mhl_sub<median(mhl_sub))]
  gene.H = names(mhl_sub)[which(mhl_sub>median(mhl_sub))]
  
  GRN.L = GRN.H = GRN
  vec = rep(0,ncol(GRN))
  id = match(gene.L,colnames(GRN))
  vec[id] = 1
  GRN.L['TP53',] = vec
  vec = rep(0,ncol(GRN))
  id = match(gene.H,colnames(GRN))
  vec[id] = 1
  GRN.H['TP53',] = vec
  
  TFA_in.L=scTFA_calculation(as.matrix(data_IntronCPM),GRN.L,zscore=F,gene_weight_method=NULL)
  TFA_ex.L=scTFA_calculation(as.matrix(data_ExonCPM),GRN.L,zscore=F,gene_weight_method=NULL)
  TFA_in.H=scTFA_calculation(as.matrix(data_IntronCPM),GRN.H,zscore=F,gene_weight_method=NULL)
  TFA_ex.H=scTFA_calculation(as.matrix(data_ExonCPM),GRN.H,zscore=F,gene_weight_method=NULL)
  
  
  
  # plot IntronTFA and ExonTFA in two replicates
  Process_data = function(TFA_in,TFA_ex){
    TFA_in_rep1 = TFA_in["TP53",grep("rep1",colnames(TFA_in))]
    TFA_in_rep2 = TFA_in["TP53",grep("rep2",colnames(TFA_in))]
    TFA_in_rep1_FC = TFA_in_rep1/as.numeric(TFA_in_rep1[1])
    TFA_in_rep2_FC = TFA_in_rep2/as.numeric(TFA_in_rep2[1])
    
    TFA_ex_rep1 = TFA_ex["TP53",grep("rep1",colnames(TFA_ex))]
    TFA_ex_rep2 = TFA_ex["TP53",grep("rep2",colnames(TFA_ex))]
    TFA_ex_rep1_FC = TFA_ex_rep1/as.numeric(TFA_ex_rep1[1])
    TFA_ex_rep2_FC = TFA_ex_rep2/as.numeric(TFA_ex_rep2[1])
    
    t_chip = c(0,1,2.5,4,5,7.5)
    chip_signals = P53_chip_data_summary$chip_data_all_FC
    
    t_final = seq(0,7.5,0.5)
    chip_LI = LinearInterpolation(t_chip,chip_signals,t_final)
    TFA_in_rep1_LI = LinearInterpolation(0:12,TFA_in_rep1,t_final)
    TFA_ex_rep1_LI = LinearInterpolation(0:12,TFA_ex_rep1,t_final)
    TFE_rep1_LI = LinearInterpolation(0:12,TF_exp_rep1,t_final)
    TFA_in_rep2_LI = LinearInterpolation(0:12,TFA_in_rep2,t_final)
    TFA_ex_rep2_LI = LinearInterpolation(0:12,TFA_ex_rep2,t_final)
    TFE_rep2_LI = LinearInterpolation(0:12,TF_exp_rep2,t_final)
    
    TFA_in_rep1_LI_FC= TFA_in_rep1_LI/TFA_in_rep1_LI[1]
    TFA_ex_rep1_LI_FC = TFA_ex_rep1_LI/TFA_ex_rep1_LI[1]
    TFE_rep1_LI_FC = TFE_rep1_LI/TFE_rep1_LI[1]
    TFA_in_rep2_LI_FC = TFA_in_rep2_LI/TFA_in_rep2_LI[1]
    TFA_ex_rep2_LI_FC = TFA_ex_rep2_LI/TFA_ex_rep2_LI[1]
    TFE_rep2_LI_FC = TFE_rep2_LI/TFE_rep2_LI[1]
    
    TFA_in_mean = (TFA_in_rep1_LI + TFA_in_rep2_LI)/2
    TFA_in_mean = TFA_in_mean/TFA_in_mean[1]
    TFA_ex_mean = (TFA_ex_rep1_LI + TFA_ex_rep2_LI)/2
    TFA_ex_mean = TFA_ex_mean/TFA_ex_mean[1]
    
    df = data.frame(chip_LI=chip_LI,TFA_in_rep1_LI_FC=TFA_in_rep1_LI_FC,TFA_ex_rep1_LI_FC=TFA_ex_rep1_LI_FC,
                    TFA_in_rep2_LI_FC=TFA_in_rep2_LI_FC,TFA_ex_rep2_LI_FC=TFA_ex_rep2_LI_FC,
                    TFA_in_mean=TFA_in_mean,TFA_ex_mean=TFA_ex_mean)
    return(df)
  }
  df1 = Process_data(TFA_in.L,TFA_ex.L)
  df2 = Process_data(TFA_in.H,TFA_ex.H)
  
  pdf(paste(result_dir,"effect_of_mRNA_halflife",".pdf",sep = ""),height=5,width=5)
  
  par(mfrow = c(1,1))
  plot(t_final,df1$TFA_in_rep1_LI_FC,col = "red",type = "b",log = "y",ylim = c(0.8,5),main = "rep1 dynamics (short mRNA half-life)",xlab = "time (h)",ylab = 'Fold change')
  points(t_final,df1$TFA_ex_rep1_LI_FC,col = "blue",type = "b")
  points(t_final,df1$chip_LI,col = "purple",type = "b")

  plot(t_final,df1$TFA_in_rep2_LI_FC,col = "red",type = "b",log = "y",ylim = c(0.8,5),main = "rep2 dynamics (short mRNA half-life)",xlab = "time (h)",ylab = 'Fold change')
  points(t_final,df1$TFA_ex_rep2_LI_FC,col = "blue",type = "b")
  points(t_final,df1$chip_LI,col = "purple",type = "b")
  
  
  plot(t_final,df2$TFA_in_rep1_LI_FC,col = "red",type = "b",log = "y",ylim = c(0.8,5),main = "rep1 dynamics (long mRNA half-life)",xlab = "time (h)",ylab = 'Fold change')
  points(t_final,df2$TFA_ex_rep1_LI_FC,col = "blue",type = "b")
  points(t_final,df2$chip_LI,col = "purple",type = "b")
  
  plot(t_final,df2$TFA_in_rep2_LI_FC,col = "red",type = "b",log = "y",ylim = c(0.8,5),main = "rep2 dynamics (long mRNA half-life)",xlab = "time (h)",ylab = 'Fold change')
  points(t_final,df2$TFA_ex_rep2_LI_FC,col = "blue",type = "b")
  points(t_final,df2$chip_LI,col = "purple",type = "b")

  cors = c(cor(df1$TFA_in_rep1_LI_FC,df1$chip_LI),cor(df1$TFA_ex_rep1_LI_FC,df1$chip_LI),
           cor(df1$TFA_in_rep2_LI_FC,df1$chip_LI),cor(df1$TFA_ex_rep2_LI_FC,df1$chip_LI),
           cor(df2$TFA_in_rep1_LI_FC,df2$chip_LI),cor(df2$TFA_ex_rep1_LI_FC,df2$chip_LI),
           cor(df2$TFA_in_rep2_LI_FC,df2$chip_LI),cor(df2$TFA_ex_rep2_LI_FC,df2$chip_LI))
  par(mfrow = c(1,1))
  all_cols = c("red","blue")
  col_ids = c(1,2,1,2,1,2,1,2)
  cols = all_cols[col_ids]
  
  names(cors) = c("rep1.L","rep1.L","rep2.L","rep2.L","rep1.H","rep1.H","rep2.H","rep2.H")
  barplot(cors,main = "correlation with ChIP signal",col = cols)
  
  # mean results
  cors = c(cor(df1$TFA_in_mean,df1$chip_LI),cor(df1$TFA_ex_mean,df1$chip_LI),
           cor(df2$TFA_in_mean,df1$chip_LI),cor(df2$TFA_ex_mean,df1$chip_LI))
  all_cols = c("red","blue")
  col_ids = c(1,2,1,2)
  cols = all_cols[col_ids]
  names(cors) = c("L","L","H","H")
  barplot(cors,main = "correlation with ChIP signal",col = cols)
  dev.off()
}


