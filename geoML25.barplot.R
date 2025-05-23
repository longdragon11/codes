#install.packages("reshape2")
#install.packages("ggpubr")
#install.packages("corrplot")


#引用包
library(reshape2)
library(ggpubr)
library(corrplot)

inputFile="CIBERSORT-Results.txt"     #免疫细胞浸润的结果文件
setwd("C:\\Users\\lexb\\Desktop\\geoML\\25.barplot")     #设置工作目录

#读取免疫细胞浸润文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

#对样品进行分组(对照组和实验组)
con=grepl("_Control", rownames(rt), ignore.case=T)
treat=grepl("_Treat", rownames(rt), ignore.case=T)
conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
data=t(rbind(conData, treatData))

#绘制柱状图
pdf(file="barplot.pdf", width=13, height=7)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data,col=col,xaxt="n",yaxt="n",ylab="Relative Percent",cex.lab=1.5)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
#在图形中标注对照组的样品,用绿色表示
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="#008B45FF")
text(a1[conNum]/2,-0.035,"Control",cex=1.8)
#在图形中标注实验组的样品,用红色表示
rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5, ytop = -0.06,col="#EE0000FF")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"Treat",cex=1.8)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1)
dev.off()

##################绘制箱线图##################
#把数据转换成ggplot2输入文件
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=cbind(as.data.frame(t(data)), Type)
data=melt(data, id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")
#绘制箱线图
group=levels(factor(data$Type))
bioCol=c("#008B45FF","#EE0000FF","#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Type",
				  xlab="",
				  ylab="Fraction",
				  legend.title="Type",
				  #notch=T, add="point",
				  width=0.8,
				  palette=bioCol)+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#输出箱线图
pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()

##################相关性图##################
treatData=treatData[,apply(treatData,2,sd)>0]
pdf(file="corHeatmap.pdf", width=12, height=12)
corrplot(corr=cor(treatData, method="spearman"),
         method = "color",        #图形的展示形式
         order = "hclust",        #免疫细胞的排序方式
         tl.col="black",          #字体颜色
         number.cex = 0.8,        #相关系数字体大小
         addCoef.col = "black",   #相关系数字体颜色
         col=colorRampPalette(c("blue", "white", "red"))(50),    #颜色的设置
         )
dev.off()

