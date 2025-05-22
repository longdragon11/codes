#install.packages("VennDiagram")


library(VennDiagram)             #引用包
diffFile="diff.txt"              #差异分析的结果文件
moduleFile="module_red.txt"      #模块基因的文件
setwd("C:\\Users\\lexb\\Desktop\\geoML\\11.venn")     #设置工作目录
geneList=list()

#读取差异分析的结果文件
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #提取差异基因的名称
geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
uniqGene=unique(geneNames)               #对差异基因取unique
geneList[["DEG"]]=uniqGene

#读取模块基因的文件
rt=read.table(moduleFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #提取模块基因的名称
geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
uniqGene=unique(geneNames)               #对模块基因取unique
geneList[["WGCNA"]]=uniqGene

#绘制venn图
venn.plot=venn.diagram(geneList,filename=NULL,fill=c("cornflowerblue", "darkorchid1"),scaled=FALSE,cat.pos=c(-1,1),cat.col = c("cornflowerblue", "darkorchid1"),cat.cex = 1.2)
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#输出交集基因的列表
interGenes=Reduce(intersect, geneList)
write.table(interGenes, file="interGenes.txt", sep="\t", quote=F, col.names=F, row.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: seqBio


