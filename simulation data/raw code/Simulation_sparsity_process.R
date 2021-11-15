# process figure for plot  # 2021/11/6
# load data
cor_sp = read.csv(file = "./cor_sp.csv",header = F)
cor_un = read.csv(file = "./cor_un.csv",header = F)
pinfo = read.csv(file = "./pinfo.csv",header = F)
cor_ratio = cor_un/cor_sp
n = sqrt(nrow(cor_sp))
mat = matrix(0,nrow = n,ncol = n)
for(i in 1:nrow(mat)){
  for(j in 1:ncol(mat)){
    id = (i-1)*n+j
    mat[i,j] = cor_ratio[id,1]
  }
}
rownames(mat) = seq(0.05,0.3,0.05)
colnames(mat) = seq(0.05,0.3,0.05)

mat = mat[n:1,]
library(pheatmap)
pheatmap(mat,cluster_rows = F,cluster_cols = F)
pdf(file ="../results/heatmap_sparsity.pdf",width = 6,height = 5)
pheatmap(mat,cluster_rows = F,cluster_cols = F)
dev.off()
