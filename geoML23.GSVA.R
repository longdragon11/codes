#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("GSEABase")
#BiocManager::install("GSVA")

#install.packages("ggpubr")


#引用包
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

gene="FIBP"      #基因的名称(根据研究的目标基因进行修改)
expFile="merge.normalize.txt"           #表达数据文件
gmtFile="c2.cp.kegg.Hs.symbols.gmt"     #基因集文件
setwd("C:\\Users\\lexb\\Desktop\\geoML\\23.GSVA")     #设置工作目录

#读取表达输入文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#去除对照组的样品
Type=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data=data[,Type=="Treat",drop=F]

#读取基因集文件
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#GSVA分析
gsvaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#对打分进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
gsvaScore=normalize(gsvaScore)
gsvaScore=gsvaScore[apply(gsvaScore,1,sd)>0.01,]

#根据目标基因的表达量对样品进行分组
lowName=colnames(data)[data[gene,]<median(data[gene,])]       #低表达组的样品
highName=colnames(data)[data[gene,]>=median(data[gene,])]     #高表达组的样品
lowScore=gsvaScore[,lowName]
highScore=gsvaScore[,highName]
data=cbind(lowScore, highScore)
conNum=ncol(lowScore)
treatNum=ncol(highScore)
Type=c(rep("Control",conNum), rep("Treat",treatNum))

#对通路进行循环, 通路差异分析
outTab=data.frame()
for(i in row.names(data)){
	test=t.test(data[i,] ~ Type)
	t=test$statistic
	pvalue=test$p.value
	if(test$estimate[2]>test$estimate[1]){t=abs(t)}else{t=-abs(t)}
	Sig=ifelse(pvalue>0.05, "Not", ifelse(t>0,"Up","Down"))
	outTab=rbind(outTab, cbind(Pathway=i, t, pvalue, Sig))
}
notSigTab=outTab[outTab$Sig=="Not",]
notSigTab=notSigTab[order(as.numeric(notSigTab$t)),]
sigTab=outTab[outTab$Sig!="Not",]
sigTab=sigTab[order(as.numeric(sigTab$t)),]
if(nrow(sigTab)>20){
	outTab=rbind(sigTab[c(1:10,((nrow(sigTab)-9):nrow(sigTab))),],notSigTab[c(1:5,((nrow(notSigTab)-4):nrow(notSigTab))),])
}else{
	outTab=rbind(sigTab,notSigTab[c(1:5,((nrow(notSigTab)-4):nrow(notSigTab))),])
}

#绘制柱状图
pdf(file="barplot.pdf", width=10.5, height=7)
outTab$t=as.numeric(outTab$t)
outTab$Sig=factor(outTab$Sig, levels=c("Down", "Not", "Up"))
gg1=ggbarplot(outTab, x="Pathway", y="t", fill = "Sig", color = "white",
		palette=c("green3","grey","red3"), sort.val = "asc", sort.by.groups = T,
		rotate=TRUE, title=gene, legend.title="Group", legend="right",
		xlab="", ylab="t value of GSVA score", x.text.angle=60)
print(gg1)
dev.off()



