#install.packages("ggplot2")
#install.packages("ggrepel")


#引用包
library(dplyr)
library(ggplot2)
library(ggrepel)

logFCfilter=0.585              #logFC过滤条件
adj.P.Val.Filter=0.05          #矫正后的p值过滤条件
diffFile="all.txt"             #所有基因差异分析的结果文件
geneFile="model.genes.txt"     #模型基因的文件
method="glmBoost+RF"           #选择机器学习的方法(需要根据热图进行修改)
setwd("C:\\Users\\lexb\\Desktop\\geoML\\18.vol")       #设置工作目录

#读取差异分析的结果文件
rt=read.table(diffFile, header=T, sep="\t", check.names=F)
row.names(rt)=rt[,1]
#定义显著性
Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")

#绘制火山图
rt = mutate(rt, Sig=Sig)
p = ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
    geom_point(aes(col=Sig))+
    scale_color_manual(values=c("green", "grey","red"))+
    labs(title = " ")+
    theme(plot.title = element_text(size=16, hjust=0.5, face = "bold"))

#在图形中标注特征基因的名称
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
geneRT=geneRT[geneRT$algorithm==method,]
sameGene=intersect(as.vector(geneRT[,1]), row.names(rt))
showData=rt[sameGene,]
p1=p+geom_label_repel(data=showData,
                    box.padding=0.2, point.padding=0.2, min.segment.length=0.1,
                    size=3, aes(label=id)) + theme_bw()
#输出图形
pdf(file="vol.pdf", width=5.25, height=4.5)
print(p1)
dev.off()

#输出最优模型基因的列表
write.table(sameGene, file="modelGene.list.txt", sep="\t", quote=F, row.names=F, col.names=F)
#输出最优模型基因的差异情况
write.table(showData, file="modelGene.diff.txt", sep="\t", quote=F, row.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


