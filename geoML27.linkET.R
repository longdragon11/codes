#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("tidyverse")

#install.packages("devtools")
#devtools::install_github("Hy4m/linkET", force = TRUE)


#引用包
library(limma)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(linkET)

expFile="merge.normalize.txt"       #表达数据文件
geneFile="modelGene.list.txt"       #基因列表文件
immFile="CIBERSORT-Results.txt"     #免疫细胞浸润的结果文件
setwd("C:\\Users\\lexb\\Desktop\\geoML\\27.linkET")     #设置工作目录

#读取表达数据文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#读取基因列表文件, 提取模型基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

#去除对照组样品
group=gsub("(.*)\\_(.*?)", "\\", colnames(data))
data=data[,group=="Treat",drop=F]
data=t(data)

#读取免疫细胞结果文件，并对数据进行整理
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
#对表达数据和免疫细胞数据取交集
sameSample=intersect(row.names(data), row.names(immune))
data=data[sameSample,,drop=F]
immune=immune[sameSample,,drop=F]
immune=immune[,apply(immune,2,sd)>0]
geneLists=list()
for(i in 1:ncol(data)){geneLists[[colnames(data)[i]]]=i}

#基因与免疫细胞相关性分析
geneCor=data.frame()
for(cell in colnames(immune)){
	if(sd(immune[,cell])==0){next}
	for(gene in colnames(data)){
		x=as.numeric(immune[,cell])
		y=as.numeric(data[,gene])
		corT=cor.test(x, y, method="spearman")
		cor=corT$estimate
		pvalue=corT$p.value
		geneCor=rbind(geneCor, cbind(spec=gene, env=cell, r=cor, p=pvalue))
	}
}
geneCor$r=as.numeric(geneCor$r)
geneCor$p=as.numeric(geneCor$p)
geneCor$pd=ifelse(geneCor$p<0.05, ifelse(geneCor$r>0, "Postive", "Negative"), "Not")
geneCor$r=abs(geneCor$r)
geneCor=geneCor %>% mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, 0.6, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", "0.4 - 0.6",">= 0.6")))

#绘制图形
qcorPlot=qcorrplot(correlate(immune, method="spearman"), type = "lower", diag = FALSE) +
  geom_square() +
  #geom_mark(sep = '\n', size = 1, sig_level = c(0.05, 0.01, 0.001), sig_thres = 0.05, color="black") +
  geom_couple(aes(colour = pd, size = rd), 
              data = geneCor, 
              curvature = nice_curvature()) +
  #设置图形的颜色和图例的名称
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu"))) +
  scale_size_manual(values = c(0.5, 1.5, 2, 3)) +
  scale_colour_manual(values = c("#1B9E77", "#CCCCCC99", "#D95F02")) +
  guides(size = guide_legend(title = "abs(Cor)",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "pvalue", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Cell-cell cor", order = 3))

#输出图形
pdf(file="linkET.pdf", width=9, height=7)
print(qcorPlot)
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio

