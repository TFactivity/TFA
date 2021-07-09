rm(list = ls())

result_dir="./"
#import GRN
human_dorothea_logic_matrix=read.csv("../../../materials/human_GRN_classAB_activation_only/dorothea_human_AB.csv",
                                     row.names = 1,check.names=FALSE)

#import MCF7 half lives data
MCF7_half_life_data=read.table("./GSE49831_MCF7_halflives.txt",header=T)
MCF7_half_life_data$MCF7_half_life_mean=apply(MCF7_half_life_data,1,mean)

P53_target=colnames(human_dorothea_logic_matrix)[human_dorothea_logic_matrix["TP53",]==1]
P53_target_with_half_life=intersect(P53_target,rownames(MCF7_half_life_data))

P53_target_half_life=data.frame(gene=P53_target_with_half_life,
                                half_life_min=MCF7_half_life_data[P53_target_with_half_life,"MCF7_half_life_mean"])
P53_target_half_life=P53_target_half_life[order(P53_target_half_life$half_life_min),]
write.csv(P53_target_half_life,paste("P53_target_half_life",".csv",sep = ""))
