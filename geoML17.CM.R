rsFile="model.classMatrix.txt"      #分类的矩阵文件
method="glmBoost+RF"                #选择机器学习的方法(需要根据热图进行修改)
setwd("C:\\Users\\lexb\\Desktop\\geoML\\17.confusion")     #设置工作目录

#读取分类的矩阵文件
riskRT=read.table(rsFile, header=T, sep="\t", check.names=F, row.names=1)
#提取数据集的ID
CohortID=gsub("(.*)\\_(.*)\\_(.*)", "\\1", rownames(riskRT))
CohortID=gsub("(.*)\\.(.*)", "\\1", CohortID)
riskRT$Cohort=CohortID

#定义混淆矩阵的函数
draw_confusion_matrix <- function(cm=null, titleName=null) {
#	layout(matrix(c(1,1,2)))
	par(mar=c(2,2,3,2))
	plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
	title(paste0('CONFUSION MATRIX (', titleName,')'), cex.main=1.5)
	
	# create the matrix 
	rect(150, 430, 240, 370, col='#3F97D0')
	text(195, 435, 'Control', cex=1.2)
	rect(250, 430, 340, 370, col='#F7AD50')
	text(295, 435, 'Treat', cex=1.2)
	text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
	text(245, 450, 'Actual', cex=1.3, font=2)
	rect(150, 305, 240, 365, col='#F7AD50')
	rect(250, 305, 340, 365, col='#3F97D0')
	text(140, 400, 'Control', cex=1.2, srt=90)
	text(140, 335, 'Treat', cex=1.2, srt=90)
	  
	# add in the cm results 
	res <- as.numeric(cm)
	text(195, 400, res[1], cex=1.6, font=2, col='white')
	text(195, 335, res[2], cex=1.6, font=2, col='white')
	text(295, 400, res[3], cex=1.6, font=2, col='white')
	text(295, 335, res[4], cex=1.6, font=2, col='white')
}

#对数据集进行循环
for(Cohort in unique(riskRT$Cohort)){
	#提取样品的分组信息(对照组和实验组)
	rt=riskRT[riskRT$Cohort==Cohort,]
	y=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(rt))
	y=ifelse(y=="Control", 0, 1)
	
	#混淆矩阵
	result_matrix=table(rt[,method], y)
	pdf(file=paste0("confusion.", Cohort, ".pdf"), width=6, height=5)
	draw_confusion_matrix(cm=result_matrix, titleName=Cohort)
	dev.off()
}


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


