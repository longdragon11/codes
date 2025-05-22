#install.packages(c("seqinr", "plyr", "openxlsx", "randomForestSRC", "glmnet", "RColorBrewer"))
#install.packages(c("ade4", "plsRcox", "superpc", "gbm", "plsRglm", "BART", "snowfall"))
#install.packages(c("caret", "mboost", "e1071", "BART", "MASS", "pROC", "xgboost"))

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("mixOmics")
#BiocManager::install("survcomp")
#BiocManager::install("ComplexHeatmap")


#引用包
library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)


#设置工作目录
setwd("C:\\Users\\lexb\\Desktop\\geoML\\15.ML")
source("refer.ML.R")

#读取训练组的数据文件
Train_data <- read.table("data.train.txt", header = T, sep = "\t", check.names=F, row.names=1, stringsAsFactors=F)
Train_expr=Train_data[,1:(ncol(Train_data)-1),drop=F]
Train_class=Train_data[,ncol(Train_data),drop=F]

#读取验证组的数据文件
Test_data <- read.table("data.test.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors = F)
Test_expr=Test_data[,1:(ncol(Test_data)-1),drop=F]
Test_class=Test_data[,ncol(Test_data),drop=F]
Test_class$Cohort=gsub("(.*)\\_(.*)\\_(.*)", "\\1", row.names(Test_class))
Test_class=Test_class[,c("Cohort", "Type")]

#获取训练组和验证组的交集基因
comgene <- intersect(colnames(Train_expr), colnames(Test_expr))
Train_expr <- as.matrix(Train_expr[,comgene])
Test_expr <- as.matrix(Test_expr[,comgene])
Train_set = scaleData(data=Train_expr, centerFlags=T, scaleFlags=T) 
names(x = split(as.data.frame(Test_expr), f = Test_class$Cohort))
Test_set = scaleData(data = Test_expr, cohort = Test_class$Cohort, centerFlags = T, scaleFlags = T)

#读取机器学习方法的文件
methodRT <- read.table("refer.methodLists.txt", header=T, sep="\t", check.names=F)
methods=methodRT$Model
methods <- gsub("-| ", "", methods)


#准备机器学习模型的参数
classVar = "Type"         #设置分类的变量名
min.selected.var = 5      #基因数目的阈值
Variable = colnames(Train_set)
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))


######################根据训练组数据构建机器学习模型######################
#根据模型组合第一种机器学习方法筛选变量
preTrain.var <- list()       #用于保存各算法筛选的变量
set.seed(seed = 123)         #设置种子
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method,              #机器学习方法
                                 Train_set = Train_set,        #训练组的基因表达数据
                                 Train_label = Train_class,    #训练组的分类数据
                                 mode = "Variable",            #选择运行模式(筛选变量)
                                 classVar = classVar)
}
preTrain.var[["simple"]] <- colnames(Train_set)

#根据模型组合第二种机器学习方法构建模型
model <- list()            #初始化模型结果列表
set.seed(seed = 123)       #设置种子
Train_set_bk = Train_set
for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method
  method <- strsplit(method, "\\+")[[1]]
  if (length(method) == 1) method <- c("simple", method)
  Variable = preTrain.var[[method[1]]]
  Train_set = Train_set_bk[, Variable]
  Train_label = Train_class
  model[[method_name]] <- RunML(method = method[2],           #机器学习方法
                                Train_set = Train_set,        #训练组的基因表达数据
                                Train_label = Train_label,    #训练组的分类数据
                                mode = "Model",               #选择运行模式(构建模型)
                                classVar = classVar)
  
  #如果某种机器学习方法筛选出的变量小于阈值，则该方法结果为空
  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}
Train_set = Train_set_bk; rm(Train_set_bk)
#保存所有机器学习模型的结果
saveRDS(model, "model.MLmodel.rds")

#构建多变量逻辑回归模型
FinalModel <- c("panML", "multiLogistic")[2]
if (FinalModel == "multiLogistic"){
  logisticmodel <- lapply(model, function(fit){    #根据逻辑回归模型计算每个样本的分类概率
    tmp <- glm(formula = Train_class[[classVar]] ~ .,
               family = "binomial", 
               data = as.data.frame(Train_set[, ExtractVar(fit)]))
    tmp$subFeature <- ExtractVar(fit)
    return(tmp)
  })
}
#保存最终以多变量逻辑回归模型
saveRDS(logisticmodel, "model.logisticmodel.rds")


#根据基因表达量计算每个样本的分类得分
model <- readRDS("model.MLmodel.rds")            #使用各个机器学习模型的线性组合函数计算得分
#model <- readRDS("model.logisticmodel.rds")     #使用多变量逻辑回归模型计算得分
methodsValid <- names(model)                     #根据特征基因数目提取有效的模型
#根据基因表达量预测样本的风险得分
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], new_data = rbind.data.frame(Train_set,Test_set))
}
riskTab=as.data.frame(t(do.call(rbind, RS_list)))
riskTab=cbind(id=row.names(riskTab), riskTab)
write.table(riskTab, "model.riskMatrix.txt", sep="\t", row.names=F, quote=F)

#根据基因表达量预测样品的分类
Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(fit = model[[method]], new_data = rbind.data.frame(Train_set,Test_set))
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
#Class_mat <- cbind.data.frame(Test_class, Class_mat[rownames(Class_mat),]) # 若要合并测试集本身的样本信息文件可运行此行
classTab=cbind(id=row.names(Class_mat), Class_mat)
write.table(classTab, "model.classMatrix.txt", sep="\t", row.names=F, quote=F)

#提取每种机器学习方法筛选到的变量
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}
fea_df <- lapply(model, function(fit){
  data.frame(ExtractVar(fit))
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, file="model.genes.txt", sep = "\t", row.names = F, col.names = T, quote = F)

#计算每个模型的AUC值
AUC_list <- list()
for (method in methodsValid){
  AUC_list[[method]] <- RunEval(fit = model[[method]],      #机器学习模型
                                Test_set = Test_set,        #验证组的表达数据
                                Test_label = Test_class,    #验证组的分类数据
                                Train_set = Train_set,      #训练组的表达数据
                                Train_label = Train_class,  #训练组的分类数据
                                Train_name = "Train",       #训练组的标签
                                cohortVar = "Cohort",       #GEO的id
                                classVar = classVar)        #分类变量
}
AUC_mat <- do.call(rbind, AUC_list)
aucTab=cbind(Method=row.names(AUC_mat), AUC_mat)
write.table(aucTab, "model.AUCmatrix.txt", sep="\t", row.names=F, quote=F)


##############################绘制AUC热图##############################
#准备图形的数据
AUC_mat <- read.table("model.AUCmatrix.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors=F)

#根据AUC的均值对机器学习模型进行排序
avg_AUC <- apply(AUC_mat, 1, mean)
avg_AUC <- sort(avg_AUC, decreasing = T)
AUC_mat <- AUC_mat[names(avg_AUC),]
#获取最优模型(训练组+测试组的AUC均值最大)
fea_sel <- fea_list[[rownames(AUC_mat)[1]]]
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3))

#设置热图注释的颜色
CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired")
names(CohortCol) <- colnames(AUC_mat)

#绘制图形
cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(Cindex_mat = AUC_mat,       #AUC矩阵
                    avg_Cindex = avg_AUC,       #AUC均值
                    CohortCol = CohortCol,      #数据集的颜色
                    barCol = "steelblue",       #右侧柱状图的颜色
                    cellwidth = cellwidth, cellheight = cellheight,    #热图每个格子的宽度和高度
                    cluster_columns = F, cluster_rows = F)      #是否对数据进行聚类

#输出热图
pdf(file="model.AUCheatmap.pdf", width=cellwidth * ncol(AUC_mat) + 6, height=cellheight * nrow(AUC_mat) * 0.45)
draw(hm, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()


######生信自学网: https://www.biowolf.cn/
######课程链接1: https://shop119322454.taobao.com
######课程链接2: https://ke.biowolf.cn
######课程链接3: https://ke.biowolf.cn/mobile
######光俊老师邮箱: seqbio@foxmail.com
######光俊老师微信: eduBio


