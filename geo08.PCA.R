#install.packages("ggplot2")
#install.packages("ggpubr")


#引用包
library(ggplot2)
library(ggpubr)
setwd("C:\\Users\\19809\\Desktop\\MR\\pca")    #设置工作目录

#定义PCA分析的函数
bioPCA=function(inputFile=null, outFile=null, titleName=null){
	#读取输入文件,提取数据
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
	data=t(rt)
	Project=gsub("(.*?)\\_.*", "\\1", rownames(data))    #获取GEO数据库研究的id
	
	#PCA分析
	data.pca=prcomp(data)
	pcaPredict=predict(data.pca)
	PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], Type=Project)

	#绘制PCA的图形
	pdf(file=outFile, width=5.5, height=4.25)
	p1=ggscatter(data=PCA, x="PC1", y="PC2", color="Type", shape="Type", 
	          ellipse=T, ellipse.type="norm", ellipse.border.remove=F, ellipse.alpha = 0.1,
	          size=2, main=titleName, legend="right")+
	          theme(plot.margin=unit(rep(1.5,4),'lines'), plot.title = element_text(hjust=0.5))
	print(p1)
	dev.off()
}

#调用函数, 绘制批次矫正前的图形
bioPCA(inputFile="merge.preNorm.txt", outFile="PCA.preNorm.pdf", titleName="Before batch correction")
#调用函数, 绘制批次矫正后的图形
bioPCA(inputFile="merge.normalize.txt", outFile="PCA.normalzie.pdf", titleName="After batch correction")


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


