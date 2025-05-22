#install.packages("reshape2")
#install.packages("ggplot2")


#引用包
library(reshape2)
library(ggplot2)

setwd("C:\\Users\\19809\\Desktop\\MR\\boxplot")     #设置工作目录

#定义箱线图的函数
bioBoxplot=function(inputFile=null, outFile=null, titleName=null){
	#读取输入文件,提取数据
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
	data=t(rt)
	Project=gsub("(.*?)\\_.*", "\\1", rownames(data))            #获取数据研究的id
	Sample=gsub("(.+)\\_(.+)\\_(.+)", "\\2", rownames(data))     #获取样品的名称
	data=cbind(as.data.frame(data), Sample, Project)
	
	#把数据转换成ggplot2输入文件
	rt1=melt(data, id.vars=c("Project", "Sample"))
	colnames(rt1)=c("Project","Sample","Gene","Expression")

	#绘制箱线图
	pdf(file=outFile, width=10, height=5)
	p=ggplot(rt1, mapping=aes(x=Sample, y=Expression))+
  		geom_boxplot(aes(fill=Project), notch=T, outlier.shape=NA)+
  		ggtitle(titleName)+ theme_bw()+ theme(panel.grid=element_blank())+ 
  		theme(axis.text.x=element_text(angle=45,vjust=0.5,hjust=0.5,size=2), plot.title=element_text(hjust = 0.5))
	print(p)
	dev.off()
}

#调用函数, 绘制批次矫正前的箱线图
bioBoxplot(inputFile="merge.preNorm.txt", outFile="boxplot.preNorm.pdf", titleName="Before batch correction")
#调用函数, 绘制批次矫正后的箱线图
bioBoxplot(inputFile="merge.normalize.txt", outFile="boxplot.normalzie.pdf", titleName="After batch correction")


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


