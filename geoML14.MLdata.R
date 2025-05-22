#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sva")


#引用包
library(limma)
library(sva)

geneFile="interGenes.txt"      #基因列表文件
setwd("C:\\Users\\lexb\\Desktop\\geoML\\14.MLdata")      #设置工作目录

#获取目录下所有"normalize.txt"结尾的文件
files=dir()
files=grep("normalize.txt$", files, value=T)
geneList=list()

#读取所有表达数据文件中的基因信息，保存到geneList
for(file in files){
    rt=read.table(file, header=T, sep="\t", check.names=F)      #读取输入文件
    geneNames=as.vector(rt[,1])      #提取基因名称
    uniqGene=unique(geneNames)       #基因取unique
    header=unlist(strsplit(file, "\\.|\\-"))
    geneList[[header[1]]]=uniqGene
}

#获取表达数据的交集基因
interGenes=Reduce(intersect, geneList)

#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
    inputFile=files[i]
    header=unlist(strsplit(inputFile, "\\.|\\-"))
    #读取输入文件，并对输入文件进行整理
    rt=read.table(inputFile, header=T, sep="\t", check.names=F)
    rt=as.matrix(rt)
    rownames(rt)=rt[,1]
    exp=rt[,2:ncol(rt)]
    dimnames=list(rownames(exp),colnames(exp))
    data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
    rt=avereps(data)
    colnames(rt)=paste0(header[1], "_", colnames(rt))

    #数据合并
    if(i==1){
    	allTab=rt[interGenes,]
    }else{
    	allTab=cbind(allTab, rt[interGenes,])
    }
    batchType=c(batchType, rep(i,ncol(rt)))
}

#提取交集基因的表达量
svaTab=ComBat(allTab, batchType, par.prior=TRUE)
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
geneTab=svaTab[intersect(row.names(svaTab), as.vector(geneRT[,1])),]
geneTab=t(geneTab)

#提取训练组和验证组的数据
train=grepl("^merge", rownames(geneTab), ignore.case=T)
trainExp=geneTab[train,,drop=F]
testExp=geneTab[!train,,drop=F]
rownames(trainExp)=gsub("merge_", "Train.", rownames(trainExp))
trainType=gsub("(.*)\\_(.*)\\_(.*)", "\\3", rownames(trainExp))
testType=gsub("(.*)\\_(.*)\\_(.*)", "\\3", rownames(testExp))
trainType=ifelse(trainType=="Control", 0, 1)
testType=ifelse(testType=="Control", 0, 1)
trainExp=cbind(trainExp, Type=trainType)
testExp=cbind(testExp, Type=testType)

#输出训练组和验证组的数据
trainOut=rbind(id=colnames(trainExp), trainExp)
write.table(trainOut, file="data.train.txt", sep="\t", quote=F, col.names=F)
testOut=rbind(id=colnames(testExp), testExp)
write.table(testOut, file="data.test.txt", sep="\t", quote=F, col.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


