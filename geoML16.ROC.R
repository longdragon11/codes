#install.packages("pROC")


#引用包
library(pROC)

rsFile="model.riskMatrix.txt"      #风险矩阵文件
method="glmBoost+RF"               #选择机器学习的方法(需要根据热图进行修改)
setwd("C:\\Users\\lexb\\Desktop\\geoML\\16.ROC")     #设置工作目录

#读取风险文件
riskRT=read.table(rsFile, header=T, sep="\t", check.names=F, row.names=1)
#获取数据集的ID
CohortID=gsub("(.*)\\_(.*)\\_(.*)", "\\1", rownames(riskRT))
CohortID=gsub("(.*)\\.(.*)", "\\1", CohortID)
riskRT$Cohort=CohortID

#对数据集进行循环
for(Cohort in unique(riskRT$Cohort)){
	#提取样品的分组信息(对照组和实验组)
	rt=riskRT[riskRT$Cohort==Cohort,]
	y=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(rt))
	y=ifelse(y=="Control", 0, 1)
	
	#绘制ROC曲线
	roc1=roc(y, as.numeric(rt[,method]))      #得到模型ROC曲线的参数
	ci1=ci.auc(roc1, method="bootstrap")      #得到ROC曲线下面积的波动范围
	ciVec=as.numeric(ci1)
	pdf(file=paste0("ROC.", Cohort, ".pdf"), width=5, height=4.75)
	plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=Cohort)
	text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
	dev.off()
}





