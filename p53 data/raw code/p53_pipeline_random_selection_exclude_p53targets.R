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

# add random TFs -- exclude TP53 targets
n_TFs = 1000
n_targets = rowSums(human_dorothea_logic_matrix)["TP53"]
TF_targets = which(human_dorothea_logic_matrix["TP53",]==1)
rand_TFs = paste("rand",1:n_TFs,sep = "_")
set.seed(0)
mat = matrix(0,nrow = n_TFs,ncol = ncol(human_dorothea_logic_matrix))
for(i in 1:n_TFs){
  rand_targets = sample(setdiff(1:ncol(human_dorothea_logic_matrix),TF_targets),n_targets)
  mat[i,rand_targets] = 1
}
rownames(mat) = rand_TFs
colnames(mat) = colnames(human_dorothea_logic_matrix)
human_dorothea_logic_matrix = rbind(human_dorothea_logic_matrix,mat)


TFA_IntronCPM=scTFA_calculation(as.matrix(data_IntronCPM),human_dorothea_logic_matrix,zscore=F,
                                                   gene_weight_method=NULL)
data_ExonCPM=data_all_combined_GRN_subset[,grep("ExonCPM",colnames(data_all_combined_GRN_subset),value = T)]
data_ExonCPM=data_ExonCPM[,-(grep("t=24 h",colnames(data_ExonCPM)))] #remove 24h data
TFA_ExonCPM=scTFA_calculation(as.matrix(data_ExonCPM),human_dorothea_logic_matrix,zscore=F,
                                                   gene_weight_method=NULL)

TFA_all_combined=cbind(TFA_IntronCPM,TFA_ExonCPM[rownames(TFA_IntronCPM),])



# import P53 ChIP-seq data
P53_chip_data=as.data.frame(read_xlsx("../raw data/41594_2017_BFnsmb3452_MOESM4_ESM.xlsx",sheet = 1))
# subset data to preserve peaks within 2kb of TSS
P53_chip_data_subset_2kb=P53_chip_data[abs(P53_chip_data$DisttoTSS)<=2000,]
P53_target=colnames(human_dorothea_logic_matrix)[human_dorothea_logic_matrix["TP53",]==1]
P53_target_in_chip=intersect(P53_target,P53_chip_data_subset_2kb$ClosestGene)

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

# plot results
if(T){
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
  p = p + theme(plot.title = element_text(size = 20))
  print(p)
  dev.off()
  
  cor_TFA_in_p53 = cor(df.mean)["chip",1]
  cor_TFA_ex_p53 = cor(df.mean)["chip",2]
  
}

# show total TFA results
TFA_to_mean = (TFA_ex_rep1_LI + TFA_ex_rep2_LI+TFA_in_rep1_LI + TFA_in_rep2_LI)/2
TFA_to_mean = TFA_to_mean/TFA_to_mean[1]
plot(t_final,TFA_in_mean,col = "red",type = "b",log = "y",ylim = c(0.8,5),main = "P53 Dynamics"
     ,xlab = "Time (h)",ylab = "Fold Change")
points(t_final,TFA_ex_mean,col = "blue",type = "b")
points(t_final,chip_LI,col = "purple",type = "b")
points(t_final,TFA_to_mean,col = "green",type = "b")


# investigate calculate random results
df.results = data.frame()
for(i in 1:1000){
  rand_TF = rand_TFs[i]
  # compare TFA with TF expression
  # expression of p53
  TF_exp_rep1 = data_ExonCPM["TP53",grep("rep1",colnames(data_ExonCPM))]
  TF_exp_rep2 = data_ExonCPM["TP53",grep("rep2",colnames(data_ExonCPM))]
  TF_exp_rep1_FC = TF_exp_rep1/as.numeric(TF_exp_rep1[1])
  TF_exp_rep2_FC = TF_exp_rep2/as.numeric(TF_exp_rep2[1])
  
  # plot IntronTFA and ExonTFA in two replicates
  TFA_in_rep1 = TFA_IntronCPM[rand_TF,grep("rep1",colnames(TFA_IntronCPM))]
  TFA_in_rep2 = TFA_IntronCPM[rand_TF,grep("rep2",colnames(TFA_IntronCPM))]
  TFA_in_rep1_FC = TFA_in_rep1/as.numeric(TFA_in_rep1[1])
  TFA_in_rep2_FC = TFA_in_rep2/as.numeric(TFA_in_rep2[1])
  
  TFA_ex_rep1 = TFA_ExonCPM[rand_TF,grep("rep1",colnames(TFA_ExonCPM))]
  TFA_ex_rep2 = TFA_ExonCPM[rand_TF,grep("rep2",colnames(TFA_ExonCPM))]
  TFA_ex_rep1_FC = TFA_ex_rep1/as.numeric(TFA_ex_rep1[1])
  TFA_ex_rep2_FC = TFA_ex_rep2/as.numeric(TFA_ex_rep2[1])
  
  # Linear Interpolation
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
  TFE_mean = (TFE_rep1_LI + TFE_rep2_LI)/2
  TFE_mean = TFE_mean/TFE_mean[1]
  
  df.mean = data.frame(TFA.in = TFA_in_mean,TFA.ex = TFA_ex_mean, TFE = TFE_mean,chip = chip_LI)
  cor(df.mean)["chip",1:3]
  df.results = rbind(df.results,cor(df.mean)["chip",1:3])
  print(i)
}
colnames(df.results) = c("TFA_in","TFA_ex","TFA_exp")
boxplot(list(df.results$TFA_in,df.results$TFA_ex),names = c("TFA_in","TFA_ex"),ylab = "correlation with ChIP signal",
        col = c("red","blue"))
pdf(file =paste(result_dir,"random_regulon_results",".pdf",sep = ""),width = 5,height = 5)

plot(density(df.results$TFA_ex),col = "blue",xlim = c(-1,1),xlab = "correlation",main ="Select regulon randomly (n=1000)")
lines(density(df.results$TFA_in),col = "red")

# abline(v = df.results$TFA_exp)
abline(v = cor_TFA_in_p53,col ="red")
abline(v = cor_TFA_ex_p53,col ="blue")
dev.off()
1 - length(which(df.results$TFA_in>cor_TFA_in_p53))/length(df.results$TFA_in)
1 - length(which(df.results$TFA_ex>cor_TFA_ex_p53))/length(df.results$TFA_ex)
id = which(df.results$TFA_in>cor_TFA_in_p53)
# compare those high correlation gene sets with dorothea genesets
high_cor_TFs = paste("rand",id,sep = "_")
GRN = human_dorothea_logic_matrix
GRN_sub = GRN[c("TP53",high_cor_TFs),]
# overlap with dorothea
num_over = NULL
for(i in 1:length(id)){
  TF = paste("rand",id[i],sep = "_")
  p53_targets = which(GRN_sub["TP53",]==1)
  rand_targets = which(GRN_sub[TF,]==1)
  over_targets = intersect(p53_targets,rand_targets)
  p = length(p53_targets)/ncol(GRN_sub)
  p*length(p53_targets)
  
  num_over = c(num_over,length(over_targets))
}
hist(num_over,main = "Overlap with dorothea GRN")
abline(v = p*length(p53_targets),col = "grey")


