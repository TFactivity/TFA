rm(list = ls())
library(pheatmap)
library(viridis)       
library(clusterProfiler)
library(org.Hs.eg.db)
library(RcisTarget)
source("../../materials//scTFA_calculation.R")
library(ggplot2)

themo_demo= theme(
  text = element_text(size=30),
  panel.background = element_rect(fill = "transparent",colour = NA), 
  panel.grid.minor = element_blank(), 
  plot.title = element_text(hjust = 0.5),
  axis.text=element_text(color='black'),
  axis.title.x = element_text(size=25),
  axis.title.y = element_text(size=25),
  axis.text.x = element_text(size=22),
  axis.text.y = element_text(size=35))

#directory of TPMcalculator results
data_dir="../raw data/analyzed20200802_2/"

#directory for storing results
result_dir="../results/"
dir.create(paste(result_dir,sep=""))

#import GRN
human_dorothea_logic_matrix=read.csv("../../materials/human_GRN_classAB_activation_only/dorothea_human_AB.csv",
                                     row.names = 1,check.names=FALSE)

#import TPMcalculator data
out_file_name_vector=list.files(data_dir)
out_file_name_vector=grep("genes.out$",out_file_name_vector,value = T)
out_file_name_vector=out_file_name_vector[order(as.numeric(sub("J","",sub("-.*$","",out_file_name_vector))))]
index=0
for(out_file_i in out_file_name_vector){
  data_i=read.table(paste(data_dir,out_file_i,sep=""),fill = T,header=T)
  data_i=na.omit(data_i)
  data_i=data_i[!duplicated(data_i$Gene_Id),]
  rownames(data_i)=data_i$Gene_Id
  
  data_i$ExonCPM=1e6*(data_i$ExonReads)/sum(data_i$Reads)
  data_i$IntronCPM=1e6*(data_i$IntronReads)/sum(data_i$Reads)
  data_i$totalCPM=1e6*(data_i$ExonReads+data_i$IntronReads)/sum(data_i$Reads)
  data_i_chosen=data_i[,c("Gene_Id","IntronLength","ExonCPM","IntronCPM","totalCPM")]
  
  #add sample name
  sample_name=sub("-LCG.*$","",out_file_i)
  colnames(data_i_chosen)[-1]=paste(colnames(data_i_chosen)[-1],sample_name,sep="_")
  
  if(index==0){
    data_all_combined=data_i_chosen
  }else{
    data_all_combined=merge(data_all_combined,data_i_chosen,by.x="Gene_Id",by.y="Gene_Id",all=T)
  }
  index=index+1
}  
rownames(data_all_combined)=data_all_combined$Gene_Id
data_all_combined[is.na(data_all_combined)]=0
# write.csv(data_all_combined,paste(result_dir,"data_all_combined",".csv"))

#choose genes in GRN
data_all_combined_GRN=data_all_combined[data_all_combined$Gene_Id %in% colnames(human_dorothea_logic_matrix),]

#separate Intron and Exon data, and calculate TFA
IntronLength_vector=apply(data_all_combined_GRN[,grep("IntronLength",colnames(data_all_combined_GRN))],
                          1,function(x) mean(x[!is.na(x)]))
data_IntronCPM=data_all_combined_GRN[IntronLength_vector>0,grep("IntronCPM",colnames(data_all_combined_GRN),value = T)]
#remove sample J15-2 due to strange profile
data_IntronCPM=data_IntronCPM[,!grepl("J15-2",colnames(data_IntronCPM))]
data_IntronCPM=data_IntronCPM[rowSums(data_IntronCPM)>0,]
TFA_IntronCPM=scTFA_calculation(as.matrix(data_IntronCPM),human_dorothea_logic_matrix,zscore=F,
                                                   gene_weight_method=NULL)

data_ExonCPM=data_all_combined_GRN[,grep("ExonCPM",colnames(data_all_combined_GRN),value = T)]
#remove sample J15-2 due to strange profile
data_ExonCPM=data_ExonCPM[,!grepl("J15-2",colnames(data_ExonCPM))]
data_ExonCPM=data_ExonCPM[rowSums(data_ExonCPM)>0,]
TFA_ExonCPM=scTFA_calculation(as.matrix(data_ExonCPM),human_dorothea_logic_matrix,zscore=F,
                                                 gene_weight_method=NULL)

write.csv(TFA_IntronCPM,paste(result_dir,"TFA_IntronCPM",".csv"))
write.csv(TFA_ExonCPM,paste(result_dir,"TFA_ExonCPM",".csv"))

TFA_all_combined=cbind(TFA_IntronCPM,TFA_ExonCPM[rownames(TFA_IntronCPM),])
write.csv(TFA_all_combined,paste(result_dir,"TFA_all_combined",".csv"))

#save TFs with totalCPM>0 for all samples
data_totalCPM=data_all_combined_GRN[,grep("totalCPM",colnames(data_all_combined_GRN),value = T)]
data_totalCPM=data_totalCPM[,!grepl("J15-2",colnames(data_totalCPM))]
data_totalCPM=data_totalCPM[rowSums(data_totalCPM)>0,]
TFA_totalCPM=scTFA_calculation(as.matrix(data_totalCPM),human_dorothea_logic_matrix,zscore=F,
                               gene_weight_method=NULL)
TFA_all_combined = cbind(TFA_IntronCPM,TFA_ExonCPM[rownames(TFA_IntronCPM),],TFA_totalCPM[rownames(TFA_IntronCPM),])

TFA_all_combined=TFA_all_combined[which(apply(data_totalCPM[rownames(TFA_all_combined),],1,min)>0),]

#average replicates
condition_vector=paste("J",c(0,3,6,9,12,15,30,45,60),"-",sep="")
method_vector=c("ExonCPM","IntronCPM","totalCPM")
index=0
for(condition_i in condition_vector){
  for(method_i in method_vector){
    TFA_all_combined_chosen=TFA_all_combined[,grep(condition_i,
                                                   grep(method_i,colnames(TFA_all_combined),value = T),
                                                   value=T),drop=F]
    
    TFA_summary_i=data.frame(average_TFA=apply(TFA_all_combined_chosen,1,mean),
                             cellnum=ncol(TFA_all_combined_chosen))
    colnames(TFA_summary_i)=paste(colnames(TFA_summary_i),condition_i,method_i,sep="_")
    rownames(TFA_summary_i)=rownames(TFA_all_combined_chosen)
    
    if(index==0){
      TFA_summary_combined=TFA_summary_i
    }else{
      TFA_summary_combined=cbind(TFA_summary_combined,TFA_summary_i)
    }
    index=index+1
  }
}
write.csv(TFA_summary_combined,paste(result_dir,"TFA_summary_combined",".csv"))

#calculate correlation between TFAs and perform GO analysis for TF modules
method_selected_vector=c("IntronCPM","ExonCPM","totalCPM")
last_TF_in_first_module_vector=c("SMAD4","NFKB2","MYC")
last_TF_in_first_module_vector=c("SMAD4","TFAP2A","MYC")
last_TF_in_first_module_vector=c("SMAD3", "RUNX1", "WT1")

options(warn = -1)
for(i_method_selected_vector in 1:length(method_selected_vector))
{
  method_selected=method_selected_vector[i_method_selected_vector]
  
  TFA_summary_combined_method_i=TFA_summary_combined[,paste("average_TFA",condition_vector,method_selected,sep="_")]
  # remove 60min result
  #TFA_summary_combined_method_i = TFA_summary_combined_method_i[,1:8]
  
  pdf(file =paste(result_dir,"Pearson corr heatmap of TFA",method_selected,"by time",".pdf"),width = 15,height = 15)
  p1=pheatmap(cor(t(TFA_summary_combined_method_i),method="pearson"),
              color = inferno(101),breaks = seq(-1,1,length.out = 101),
              main=paste(method_selected,"(pearson correlation)"))
  graphics.off()
  
  TF_name_in_heatmap=rownames(TFA_summary_combined_method_i)[p1$tree_row$order]
  
  pdf(file =paste(result_dir,"TFA_summary",method_selected,".pdf"),width = 3,height = 15)
  pheatmap(TFA_summary_combined_method_i[TF_name_in_heatmap,],scale="row",
           cluster_rows = F,cluster_cols = F,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(101),
           main=paste("TFA_summary",method_selected))
  graphics.off()
  
  write.csv(TFA_summary_combined_method_i[TF_name_in_heatmap,],file= paste(result_dir,"TFA_summary",method_selected,".csv"),
            quote = F,row.names = T)
  pdf(file =paste(result_dir,"Two_genes_dynamics",method_selected,".pdf"),width = 8,height = 4)
  genes = c("IRF9","GATA3")
  par(mfrow = c(1,2))
  TFA_summary_combined_method_i[genes,]
  for(gene in genes){
    time_points = c(0,3,6,9,12,15,30,45,60)
    plot(time_points,TFA_summary_combined_method_i[gene,],main = paste(method_selected,gene),type ="b",
         xlab = "t(min)",ylab = "TFA (CPM)")
    
  }
  graphics.off()
  

  
  last_TF_in_first_module=last_TF_in_first_module_vector[i_method_selected_vector]
  TF_module_list=list(module1=TF_name_in_heatmap[1:which(TF_name_in_heatmap==last_TF_in_first_module)],
                      module2=TF_name_in_heatmap[(which(TF_name_in_heatmap==last_TF_in_first_module)+1):length(TF_name_in_heatmap)])
  
  #retrieve TF list from RcisTarget package
  data(motifAnnotations_hgnc)
  human_TF_list = unique(motifAnnotations_hgnc$TF)
  
  print(paste(method_selected))
  print(TF_module_list)
  if(F){
    for(i_clust in 1:length(TF_module_list)){
      TF_clust_i=TF_module_list[[i_clust]]
      
      # tmp = bitr(TF_clust_i, fromType="ALIAS", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      
      # SYMBOL -> ENSEMBLE -> ENTREZID
      # TFs = rownames(TFA_summary_combined_method_i)
      # n = length(TFs)
      # df.tmp = NULL
      # tmp1 = bitr(TFs, fromType="ALIAS", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      # 
      # tmp2 = bitr(TFs, fromType="ALIAS", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
      # 
      # tmp3 = bitr(TFs, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      # tmp4 = bitr(TFs, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db")
      # df.tmp = data.frame(from = c("ALIAS","ALIAS","SYMBOL","SYMBOL"), num_from = c(n,n,n,n),
      #                     to = c("ENTREZID","ENSEMBL","ENTREZID","ENSEMBL"),
      #                     num_to= c(nrow(tmp1),nrow(tmp2),nrow(tmp3),nrow(tmp4)))
      
      bitr_i=bitr(TF_clust_i, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")[,"ENTREZID"]
      bg_bitr=bitr(human_TF_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")[,"ENTREZID"]
      bitr_GO <- enrichGO(
        gene = bitr_i, keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
        pvalueCutoff = 1, qvalueCutoff  = 1, readable = F,
        universe=bg_bitr)
      # bitr_i = TF_clust_i
      # bg_bitr = human_TF_list
      # bitr_GO <- enrichGO(
      #   gene = bitr_i, keyType = "SYMBOL", OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
      #   pvalueCutoff = 1, qvalueCutoff  = 1, readable = F,
      #   universe=bg_bitr)
      write.csv(bitr_GO@result,paste(result_dir,"bitr_GO@result cluster",i_clust,method_selected,".csv"))
      print(paste(method_selected,'module',i_clust))
      print(bitr_GO@result[1:5,1:6])
    } 
  }
  
}

#retrieve GO list in human cells from GO.db
library(GO.db)
GOTERM_all=as.list(GOTERM)
#select GO of BP
GOTERM_all_is_BP=unlist(lapply(GOTERM_all,function(x) Ontology(x)=="BP"))
GOTERM_all_BP=GOTERM_all[GOTERM_all_is_BP]
GOTERM_all_BP=names(GOTERM_all_BP)

#retrieve gene list of each GO based on biomaRt package
library(biomaRt)
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

TF_vector=rownames(TFA_summary_combined)
TF_vector=TF_vector[!duplicated(TF_vector)]
gene_GO_all=getBM(c("hgnc_symbol","go_id"), filters = c("go","hgnc_symbol"), values = list(GOTERM_all_BP,TF_vector), mart = ensembl)

TF_num_thre=5 #preserve GOs with no less than 5 TFs with TFA data
GO_vector_chosen=gene_GO_all$go_id[!duplicated(gene_GO_all$go_id)]
GO_gene_list=list()
for(GO_vector_chosen_i in GO_vector_chosen){
  gene_i=gene_GO_all[gene_GO_all$go_id==GO_vector_chosen_i,"hgnc_symbol"]
  gene_i=gene_i[!duplicated(gene_i)]
  if(length(gene_i)>=TF_num_thre){ 
    GO_gene_list=c(GO_gene_list,list(gene_i))
    names(GO_gene_list)[length(GO_gene_list)]=GO_vector_chosen_i
  }
}
GO_gene_list=GO_gene_list[!duplicated(GO_gene_list)]

cor_computation_within_GO=function(TFA_matrix,TF_list_within_GO,correlation_method="pearson"){
  #compute correlation between TFs within one GO
  cor_within_GO=as.data.frame(t(combn(TF_list_within_GO,2)))
  cor_within_GO$cor=apply(cor_within_GO,1,function(x) cor(TFA_matrix[x[1],],TFA_matrix[x[2],],method=correlation_method))
  cor_within_GO_mean=mean(cor_within_GO$cor)
  
  return(list(cor_within_GO=cor_within_GO,
              cor_within_GO_mean=cor_within_GO_mean))
}

cor_GO_matrix=data.frame(GO=NA,TF_num_in_GO=NA,
                         cor_within_GO_mean_intron=NA,
                         cor_within_GO_mean_exon=NA)
for(i_GO_gene_list in 1:length(GO_gene_list)){
  cor_i_intron=cor_computation_within_GO(as.matrix(TFA_summary_combined[,paste("average_TFA",condition_vector,"IntronCPM",sep="_")]),
                                         GO_gene_list[[i_GO_gene_list]],correlation_method="pearson")
  cor_i_exon=cor_computation_within_GO(as.matrix(TFA_summary_combined[,paste("average_TFA",condition_vector,"ExonCPM",sep="_")]),
                                       GO_gene_list[[i_GO_gene_list]],correlation_method="pearson")
  cor_GO_matrix=rbind(cor_GO_matrix,c(names(GO_gene_list)[i_GO_gene_list],length(GO_gene_list[[i_GO_gene_list]]),
                                      cor_i_intron[["cor_within_GO_mean"]],cor_i_exon[["cor_within_GO_mean"]]))
}
cor_GO_matrix=cor_GO_matrix[-1,]
cor_GO_matrix$cor_within_GO_mean_intron=as.numeric(cor_GO_matrix$cor_within_GO_mean_intron)
cor_GO_matrix$cor_within_GO_mean_exon=as.numeric(cor_GO_matrix$cor_within_GO_mean_exon)
cor_GO_matrix=cor_GO_matrix[order(cor_GO_matrix$cor_within_GO_mean_intron,decreasing = T),]

#add GO name
cor_GO_matrix_appended=cbind(cor_GO_matrix,NA,NA)
colnames(cor_GO_matrix_appended)=c(colnames(cor_GO_matrix),"Term","Definition")
for(i_cor_GO_matrix_appended in 1:nrow(cor_GO_matrix_appended)){
  cor_GO_matrix_appended[i_cor_GO_matrix_appended,c("Term","Definition")]=
    c(Term(GOTERM_all[[cor_GO_matrix_appended[i_cor_GO_matrix_appended,"GO"]]]),
      Definition(GOTERM_all[[cor_GO_matrix_appended[i_cor_GO_matrix_appended,"GO"]]]))
}

write.csv(cor_GO_matrix_appended,paste(result_dir,"cor_GO_matrix_appended",".csv"))

p1<-ggplot(data=cor_GO_matrix_appended,mapping=aes(x=cor_within_GO_mean_exon,y=cor_within_GO_mean_intron))+
  geom_point(size=5)+
  labs(x = "cor_within_GO_mean_exon", y = "cor_within_GO_mean_intron",
       title = paste(nrow(cor_GO_matrix_appended),"GO","\n",
                     "p value for paired t test (y>x)",
                     format(t.test(cor_GO_matrix_appended$cor_within_GO_mean_intron,cor_GO_matrix_appended$cor_within_GO_mean_exon,alternative = "greater",paired = T)$p.value,
                            digits = 3)))+
  geom_abline(size=1,colour="red")+
  scale_x_continuous(limits=c(-0.25,1),breaks=seq(-0.25,1,0.25))+
  scale_y_continuous(limits=c(-0.25,1),breaks=seq(-0.25,1,0.25))+
  theme_bw()+themo_demo

pdf(file =paste(result_dir,"cor_GO_matrix_appended",".pdf"),width = 10,height =10)
print(p1)
dev.off()

# check T cell activation genes in the diagram 
if(T){
  plot(cor_GO_matrix_appended$cor_within_GO_mean_exon,cor_GO_matrix_appended$cor_within_GO_mean_intron, xlab = "cor.exonTFA",
       ylab = "cor.intronTFA",xlim = c(-0.25,1),ylim =  c(-0.25,1))
  abline(a = 0,b = 1,col= "green")
  GOs = c("GO:0045087","GO:0006954","GO:0071356","GO:0034097","GO:0032481","GO:0050728","GO:0019221","GO:0071345")
  GOs = c("GO:0045087","GO:0006954")
  ids = which(cor_GO_matrix_appended$GO %in% GOs)
  cor_GO_matrix_appended$Term[ids]
  points(cor_GO_matrix_appended$cor_within_GO_mean_exon[ids],cor_GO_matrix_appended$cor_within_GO_mean_intron[ids],
         col = "red",pch = 16)
}


# show TF expression
if(T){
  dim(data_ExonCPM)
  dim(TFA_IntronCPM)
  dim(TFA_ExonCPM)
  over_genes = intersect(rownames(TFA_IntronCPM),rownames(TFA_ExonCPM))
  TFE = data_ExonCPM[match(rownames(TFA_IntronCPM),rownames(data_ExonCPM)),]
  TFE = na.omit(TFE)
  TFs = rownames(TFE)
  TFE = as.matrix(TFE)
  
  TFE.m = apply(TFE,1,mean)
  table(TFE.m>1)
  
  # TFE dynamics
  time_points = c(0,0,3,3,6,6,9,9,12,12,15,30,30,45,45,60,60)
  par(mfrow = c(3,3))
  for(i in 1:9){
    TF = TFs[i]
    plot(time_points,TFE[i,],main = TF,xlab = "t (min)", ylab = "TF expression (CPM)")
    
  }
  i = 2
  # define change of TF
  range_TFE = apply(TFE, 1, function(x){return(max(x)-min(x))})
  length(which(range_TFE>1))
  TFE_sub = TFE[which(range_TFE>1),]
  FC_TFE = apply(TFE_sub, 1, function(x){return(max(x)/min(x))})
  sort(FC_TFE,decreasing = T)
  hist(log2(FC_TFE),breaks = 100)
  boxplot(log2(FC_TFE))
  
  ids = match(names(sort(FC_TFE,decreasing = T)[1:9]),TFs)
  n = 1
  ids = match(names(sort(FC_TFE,decreasing = T)[(1:9)+9*n]),TFs)
  
  par(mfrow = c(3,3))
  for(i in ids){
    TF = TFs[i]
    plot(time_points,TFE[i,],main = TF,xlab = "t (min)", ylab = "TF expression (CPM)")
  }
  
  # compare expression data in 0 - 15 min with 30-60min
  df.comp = data.frame()
  for(i in 1:length(TFs)){
    TF = TFs[i]
    TFE1 = TFE[TF,1:11]
    TFE2 = TFE[TF,12:17]
    pvalue = wilcox.test(TFE1,TFE2)$p.value
    FC = median(TFE2)/median(TFE1)
      
    df.comp = rbind(df.comp,c(median(TFE1),median(TFE2),pvalue,FC))
  }
  colnames(df.comp) = c("median_exp1","median_exp2","pvalue","FC")
  rownames(df.comp) = TFs
  plot(log2(df.comp$FC),-log2(df.comp$pvalue),main = "change of TF expression",xlab = "log2(Fold change)",ylab = "-log2(pvalue)")
  abline(h = -log2(0.05),col = "grey")
  abline(v = 0,col = "grey")
  
  table(df.comp$pvalue<0.05)
  sort( df.comp$FC[which(df.comp$pvalue<0.05)])  
  
  
  # correlation of TFE with TFA
  df.cor = data.frame()
  for(i in 1:length(TFs)){
    TF = TFs[i]
    cor.ex = cor(TFE[TF,],TFA_ExonCPM[TF,])
    cor.in = cor(TFE[TF,],TFA_IntronCPM[TF,])
    
    df.cor = rbind(df.cor,c(cor.ex,cor.in))
  }
  colnames(df.cor) = c("cor.ex","cor.in")
  plot(df.cor$cor.ex,df.cor$cor.in)
}

# show modules of intron+exon results
if(T){
  total_genes = union(rownames(data_ExonCPM),rownames(data_IntronCPM))
  data_totalCPM = matrix(0,nrow = length(total_genes),ncol = ncol(data_ExonCPM))
  rownames(data_totalCPM) = total_genes
  colnames(data_totalCPM) = colnames(data_ExonCPM)
  id1 = match(rownames(data_ExonCPM),total_genes)
  for(i in 1:length(id1)){
    for(j in 1:ncol(data_totalCPM)){
      data_totalCPM[id1[i],j] =  data_totalCPM[id1[i],j]+data_ExonCPM[i,j]
    }
  }
  id2 = match(rownames(data_IntronCPM),total_genes)
  for(i in 1:length(id1)){
    for(j in 1:ncol(data_totalCPM)){
      data_totalCPM[id2[i],j] =  data_totalCPM[id2[i],j]+data_ExonCPM[i,j]
    }
  }
  
  
}

# compare with ATAC-seq data
if(T){
  # load ATAC_seq data
  if(T){
    load(file = "./ATACseq_count.RData")
    mat = mat_new
    # normalize
    for(i in 1:ncol(mat)){
      mat[,i] = mat[,i]/colSums(mat)[i]*1000000
    }
    
    # calculate activity of each TF
    #import GRN
    GRN=read.csv("../../materials/human_GRN_classAB_activation_only/dorothea_human_AB.csv",
                 row.names = 1,check.names=FALSE)
    # calculate each TF activity
    dim(GRN)
    TFA = NULL
    for(i in 1:nrow(GRN)){
      targets = colnames(GRN)[which(GRN[i,]==1)]
      over_targets = intersect(targets,rownames(mat))
      mat.sub = mat[over_targets,]
      if(length(over_targets)==1){
        TFA = rbind(TFA,mat.sub)
      }else{
        TFA = rbind(TFA,apply(mat.sub,2,mean))
      }
      print(i/nrow(GRN))
    }
    rownames(TFA) = rownames(GRN)
    TFA_mean = data.frame(NC_m = apply(TFA[,1:2],1,mean),PIC_m = apply(TFA[,3:4],1,mean))
    
  }
  TFA_in = TFA_all_combined[,c(1,2,16,17)]
  TFA_in = TFA_all_combined[,c(1,2,12,13)]
  TFA_in.m = data.frame(in_0=apply(TFA_in[,1:2],1,mean),in_60=apply(TFA_in[,3:4],1,mean))
  
  
  TFA_mean = TFA_mean[rownames(TFA_in),]
  # compare fold change
  df.comp = data.frame(atac_fold = TFA_mean[,2]/TFA_mean[,1],in_fold = TFA_in.m[,2]/TFA_in.m[,1])
  pdf(file =paste(result_dir,"compare_with_ATAC-seq",".pdf"),width = 5,height =5)
  
  plot(df.comp$atac_fold,df.comp$in_fold,log = "xy",xlab = "ATAC-seq fold change",ylab = "intron fold change")
  abline(h = 1,col = "grey",lty = 2)
  abline(v = 1,col = "grey",lty = 2)
  cor.test(df.comp$atac_fold,df.comp$in_fold)
  dev.off()
  
  up1 = which(df.comp$atac_fold>1)
  do1 = which(df.comp$atac_fold<1)
  up2 = which(df.comp$in_fold>1)
  do2 = which(df.comp$in_fold<1)
  length(intersect(up1,up2))
  length(intersect(do1,up2))
  length(intersect(do1,do2))
  length(intersect(up1,do2))
  
}
