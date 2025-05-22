#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")


inputFile="merge.normalize.txt"      #输入文件
setwd("C:\\Users\\lexb\\Desktop\\geoML\\24.CIBERSORT")      #设置工作目录
source("geoML24.CIBERSORT.R")       #引用包

#免疫细胞浸润分析
outTab=CIBERSORT("ref.txt", inputFile, perm=1000)

#对免疫浸润结果过滤，并且保存免疫细胞浸润结果
outTab=outTab[outTab[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


