rm(list = ls())
library(ggplot2)

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
corr_data=as.data.frame(t(readxl::read_xls(paste(result_dir,"cor_mat.xls",sep=""),col_names = F)))

data = as.matrix(corr_data[2:nrow(corr_data),1:ncol(corr_data)])
data = apply(data,2,as.numeric)
colnames(data) = as.character(as.matrix(corr_data[1,]))
corr_data = as.data.frame(data)


##  attention!!
corr_data_reformat=data.frame(corr=as.numeric(c(corr_data$unspliced,corr_data$splice)),  
                              index=rep(c("unspliced TFA","spliced TFA"),each=nrow(corr_data)))
write.csv(corr_data_reformat,paste(result_dir,"corr_data_reformat",".csv",sep = ""))

p1<-ggplot(data=corr_data_reformat,mapping=aes(x=index,y=corr,colour=index,fill=index,group=index))+
  geom_jitter(size=1)+
  labs(x = "index", y = "Pearson correlation",
       title = paste(nrow(corr_data_reformat)/length(table(corr_data_reformat$index)),"cells"))+
  scale_y_continuous(limits=c(0,1))+
  theme_bw()+themo_demo

pdf(file =paste(result_dir,"corr_data_reformat.pdf",sep = ""),width = 8,height = 8)
print(p1)
dev.off()

