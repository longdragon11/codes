#install.packages("glmnet")
#install.packages("pROC")


#引用包
library(glmnet)
library(pROC)

expFile="merge.normalize.txt"      #表达数据文件
geneFile="modelGene.list.txt"      #基因列表文件
setwd("C:\\Users\\lexb\\Desktop\\geoML\\20.geneROC")    #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

#提取样品的分组信息(对照组和实验组)
y=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(rt))
y=ifelse(y=="Control", 0, 1)

#读取基因列表文件
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

#定义图形的颜色
bioCol=rainbow(nrow(geneRT), s=0.9, v=0.9)

#对基因进行循环，绘制ROC曲线
aucText=c()
k=0
for(x in as.vector(geneRT[,1])){
	k=k+1
	#绘制ROC曲线
	roc1=roc(y, as.numeric(rt[x,]))     #得到ROC曲线的参数
	if(k==1){
		pdf(file="ROC.genes.pdf", width=5, height=4.5)
		plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3)
		aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
	}else{
		plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3, add=TRUE)
		aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
	}
}
#绘制图例，得到ROC曲线下的面积
legend("bottomright", aucText, lwd=3, bty="n", cex=0.8, col=bioCol[1:(ncol(rt)-1)])
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


