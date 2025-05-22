#    install.packages("BiocManager")
#BiocManager::install("limma")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("sva")


#引用包
library(limma)
library(sva)
setwd("C:\\Users\\19809\\Desktop\\MR\\数据合并")      #设置工作目录

#获取目录下所有"normalize.txt"结尾的文件
files=dir()
files=grep("normalize.txt$", files, value=T)
geneList=list()

#读取所有txt文件中的基因信息，保存到geneList
for(file in files){
	if(file=="merge.preNorm.txt"){next}
	if(file=="merge.normalize.txt"){next}
    rt=read.table(file, header=T, sep="\t", check.names=F)      #读取输入文件
    geneNames=as.vector(rt[,1])      #提取基因名称
    uniqGene=unique(geneNames)       #基因取unique
    header=unlist(strsplit(file, "\\.|\\-"))
    geneList[[header[1]]]=uniqGene
}

#获取交集基因
interGenes=Reduce(intersect, geneList)

#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
    inputFile=files[i]
	if(inputFile=="merge.preNorm.txt"){next}
	if(inputFile=="merge.normalize.txt"){next}
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

#输出合并后的表达数据
outTab=rbind(geneNames=colnames(allTab), allTab)
write.table(outTab, file="merge.preNorm.txt", sep="\t", quote=F, col.names=F)

#对合并后数据进行批次矫正，输出批次矫正后的表达数据
outTab=ComBat(allTab, batchType, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge.normalize.txt", sep="\t", quote=F, col.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


