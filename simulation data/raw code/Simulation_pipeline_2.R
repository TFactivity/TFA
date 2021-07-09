rm(list = ls())
library(ggplot2)
library("easyGgplot2")


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

result_dir="../results/"

#import matlab results
corr_data=as.data.frame(t(readxl::read_xls(paste(result_dir,"corr_combined.xls",sep=""),col_names = F)))
data = as.matrix(corr_data[2:nrow(corr_data),2:ncol(corr_data)])
data = apply(data,2,as.numeric)
rownames(data) = as.character(corr_data[-1,1])
colnames(data) = as.character(as.matrix(corr_data[1,-1]))
corr_data = as.data.frame(data)
corr_data$name = rownames(corr_data)

corr_data_reformat=data.frame(target_name=factor(c(as.character(corr_data$name),as.character(corr_data$name)),
                                                 levels=as.character(corr_data$name)),
                              target_half_life=as.numeric(c(corr_data$`half life`,corr_data$`half life`)),
                              corr=as.numeric(c(corr_data$corr_TF_unspliced_mean,corr_data$corr_TF_spliced_mean)),
                              corr_se=as.numeric(c(corr_data$corr_TF_unspliced_se,corr_data$corr_TF_spliced_se)),
                              index=rep(c("corr_TF_unspliced","corr_TF_spliced"),each=nrow(corr_data)))
write.csv(corr_data_reformat,paste(result_dir,"corr_data_reformat.csv",sep = ""))  

p1<-ggplot(data=corr_data_reformat,mapping=aes(x=target_name,y=target_half_life))+
  geom_bar(stat="identity",position='dodge',width=0.5)+
  labs(x = "target_name", y = "target_half_life (min)",
       title = paste(nrow(corr_data),"targets"))+
  scale_x_discrete(labels=c())+
  theme_bw()+themo_demo

p2<-ggplot(data=corr_data_reformat,mapping=aes(x=target_name,y=corr,color=index,group=index,fill=index))+
  geom_bar(stat="identity",position=position_dodge(),size=0)+
  geom_errorbar(aes(ymin=corr-corr_se, ymax=corr+corr_se),
               position=position_dodge(.9),color="black",size=.8,width=0)+
  labs(x = "target_name", y = "Correlation",
       title = paste(nrow(corr_data),"targets,","error bar: se"))+
  scale_x_discrete(labels=c())+
  scale_y_continuous(limits=c(0,0.8),breaks=seq(0,0.8,0.2))+
  theme_bw()+themo_demo

pdf(file =paste(result_dir,"correlation",".pdf",sep = ""),width = 16,height = 20)
ggplot2.multiplot(p1,p2,cols=1)
dev.off()

