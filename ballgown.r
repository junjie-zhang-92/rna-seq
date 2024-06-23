rm(list=ls())   #清空列表
setwd("E:/hisat2/ballgown")  #设置路径
library(ballgown)
library(genefilter)
library(dplyr)
library(reshape)
library(ggplot2)

pheno_data = read.csv("SRR.csv")  #导入表现数据，或read.table("file.txt")
pheno_data
bg = ballgown(samples = as.vector(pheno_data$path),pData = pheno_data) #ballgown读取数据，path为样品路径
bg

bg_filt = subset(bg,"rowVars(texpr(bg))>1",genomesubset=TRUE) #过滤表达量为0的
bg_filt
bg_table = texpr(bg_filt, 'all') #载入样本数据
bg_gene_names = unique(bg_table[, 9:10]) #提取样本名字

results_transcripts=stattest(bg_filt,feature="transcript",covariate="type",getFC=TRUE,meas="FPKM") #差异分析
results_genes = stattest(bg_filt, feature="gene", covariate="type", getFC=TRUE, meas="FPKM")

results_transcripts=data.frame(geneNames=ballgown::geneNames(bg_filt),geneIDs=geneIDs(bg_filt),results_transcripts)  #转录本添加基因名字
results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id")) 

results_transcripts=arrange(results_transcripts,pval) #根据pval进行排序
results_genes=arrange(results_genes,pval)

write.csv(results_transcripts, "chrX_transcript_results.csv",row.names=FALSE) #保存数据到本地文件
write.csv(results_genes, "chrX_gene_results.csv",row.names=FALSE)

results_transcripts_0.05=subset(results_transcripts,results_transcripts$qval<0.05)  #筛选qval<0.05的转录本
results_genes_0.05=subset(results_genes,results_genes$qval<0.05)
write.csv(results_transcripts_0.05,file="chrX_transcript_0.05.csv",row.names=FALSE)  #保存数据
write.csv(results_genes_0.05,file="chrX_genes_0.05.csv",row.names=FALSE)

tropical= c('darkorange', 'dodgerblue','hotpink', 'limegreen', 'yellow','red')  #加载颜色
palette(tropical)
options(repr.plot.width = 9, repr.plot.height = 6)  #打开画板
fpkm = texpr(bg,meas="FPKM")  #
fpkm = log2(fpkm+1)
boxplot(fpkm,col=as.numeric(pheno_data$sex),las=2,ylab='log2(FPKM+1)')  #画箱型图
dev.off()

fpkm_L <- melt(fpkm)  #融合fpkm
tail(fpkm_L)
colnames(fpkm_L)=c('n','sample','value')  #提取fpkm_l文件中的两列，并加一列数字
head(fpkm_L)
group_list=as.character(pheno_data$type)  #将type单作一个文件夹
group_list
fpkm_L$group=rep(group_list,each=nrow(fpkm))  #加一列group
tail(fpkm_L)
options(repr.plot.width = 9, repr.plot.height = 5)  #打开画板
p=ggplot(fpkm_L,aes(x=sample,y=value,fill=group))+
  geom_boxplot() + 
  theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1,size = 6))
p  #画箱型图
ggplot(results_transcripts_0.05,aes(log(fc),-log(pval)))+geom_point()  #画火山散点图
