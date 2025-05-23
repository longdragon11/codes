#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")
#install.packages("PerformanceAnalytics")


#引用包
library(limma)
library(reshape2)
library(ggpubr)
library(PerformanceAnalytics)

expFile="merge.normalize.txt"      #表达数据文件
geneFile="modelGene.list.txt"      #基因列表文件
setwd("C:\\Users\\lexb\\Desktop\\geoML\\19.boxplot")     #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)

#读取基因列表文件, 提取模型基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(data))
data=t(data[sameGene,])

#获取样品的分组信息(对照组和实验组)
Type=gsub("(.*)\\_(.*?)", "\\2", row.names(data))
treatData=data[Type=="Treat",]
rt=cbind(as.data.frame(data), Type)

#将数据转换为箱线图的输入文件
data=melt(rt, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#绘制箱线图
p=ggboxplot(data, x="Gene", y="Expression", fill = "Type",
	     xlab="",
	     ylab="Gene expression",
	     legend.title="Type",
	     palette = c("#0088FF", "#FF5555"), width=0.75)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#输出图形
pdf(file="boxplot.pdf", width=6, height=4.5)
print(p1)
dev.off()

#绘制相关性图形(PerformanceAnalytics包)
pdf(file="cor.pdf", width=7, height=6.5)
chart.Correlation(treatData, histogram=TRUE, pch=19, method="pearson")
dev.off()




