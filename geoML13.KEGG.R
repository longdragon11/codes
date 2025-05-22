#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05      #p值过滤条件
adjPvalFilter=1        #矫正后的p值过滤条件

#定义图形的颜色
colorSel="p.adjust"
if(adjPvalFilter>0.05){
	colorSel="pvalue"
}

setwd("C:\\Users\\lexb\\Desktop\\geoML\\13.KEGG")      #设置工作目录
rt=read.table("interGenes.txt", header=F, sep="\t", check.names=F)     #读取交集基因的列表文件

#提取基因的名称, 将基因名称转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt, entrezIDs)
colnames(rt)[1]="gene"
rt=rt[rt[,"entrezIDs"]!="NA",]      #去除基因id为NA的基因
gene=rt$entrezID
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
kk@result$Description=gsub(" - Homo sapiens \\(human\\)", "", kk@result$Description)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$p.adjust<adjPvalFilter),]
#输出显著富集的结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#设置展示通路的数目
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#柱状图
pdf(file="barplot.pdf", width=8.5, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=100, color=colorSel)
dev.off()

#气泡图
pdf(file="bubble.pdf", width=8.5, height=7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=100, color=colorSel)
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


