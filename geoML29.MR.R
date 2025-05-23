#install.packages("R.utils")

#devtools::install_github("mrcieu/gwasglue",force = TRUE)
#BiocManager::install("VariantAnnotation")

#install.packages("devtools")
#devtools::install_github("mrcieu/gwasglue", force = TRUE)

#install.packages("remotes")
#remotes::install_github("MRCIEU/TwoSampleMR")


#引用包
library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)

exposureFile="gene.exposure.csv"      #暴露数据文件(选择基因或者蛋白数据)
geneFile="modelGene.diff.txt"         #基因列表文件
outcomeName="Alzheimer disease"       #设置图形中展示疾病的名称(需修改)
setwd("C:\\Users\\lexb\\Desktop\\geoML\\29.MR")      #设置工作目录

#读取暴露数据
rt=read_exposure_data(filename=exposureFile,
                                sep = ",",
                                snp_col = "SNP",
                                beta_col = "beta.exposure",
                                se_col = "se.exposure",
                                pval_col = "pval.exposure",
                                effect_allele_col="effect_allele.exposure",
                                other_allele_col = "other_allele.exposure",
                                eaf_col = "eaf.exposure",
                                phenotype_col = "exposure",
                                id_col = "id.exposure",
                                samplesize_col = "samplesize.exposure",
                                chr_col="chr.exposure", pos_col = "pos.exposure",
                                clump=FALSE)

#对结局数据进行循环
files=dir()                           #获取目录下所有文件
files=grep(".gz$",files,value=T)      #提取.gz结尾的文件
for(outcomeFile in files){
	#读取结局数据文件	
	outcomeID=gsub(".gz", "", outcomeFile)
	outcomeData=read_outcome_data(snps=rt$SNP,
	                 filename=outcomeFile, sep = "\t",
	                 snp_col = "rsids",
	                 beta_col = "beta",
	                 se_col = "sebeta",
	                 effect_allele_col = "alt",
	                 other_allele_col = "ref",
	                 pval_col = "pval",
	                 eaf_col = "af_alt")
	write.csv(outcomeData, file="outcome.csv", row.names=F)
	
	#对特征基因进行循环
	geneRT=read.table(geneFile, header=T, sep="\t", check.names=F, row.names=1)
	sameGene=intersect(row.names(geneRT), unique(rt$exposure))
	outTab=data.frame()
	for(expoName in sameGene){
		#提取这个基因的暴露数据
		i=paste0(expoName, "__", outcomeID)
		singleExposureFile=paste0(i, ".exposure.csv")
		exposure_set=rt[rt$exposure==expoName,]
		if(nrow(exposure_set)>=3){
			write.csv(exposure_set, file=singleExposureFile, row.names=F)
			#读取这个基因的暴露数据
			exposure_dat=read_exposure_data(filename=singleExposureFile,
		                                sep = ",",
		                                snp_col = "SNP",
		                                beta_col = "beta.exposure",
		                                se_col = "se.exposure",
		                                pval_col = "pval.exposure",
		                                effect_allele_col="effect_allele.exposure",
		                                other_allele_col = "other_allele.exposure",
		                                eaf_col = "eaf.exposure",
		                                phenotype_col = "exposure",
		                                samplesize_col = "samplesize.exposure",
		                                chr_col="chr.exposure", pos_col = "pos.exposure",
		                                clump=FALSE)
		    file.remove(singleExposureFile)
			
			#读取结局数据
			outcome_data=read_outcome_data(snps=exposure_dat$SNP,
				             filename="outcome.csv", sep = ",",
				             snp_col = "SNP",
				             beta_col = "beta.outcome",
				             se_col = "se.outcome",
				             effect_allele_col = "effect_allele.outcome",
				             other_allele_col = "other_allele.outcome",
				             pval_col = "pval.outcome",
				             eaf_col = "eaf.outcome")
			
			#将暴露数据和结局数据合并
			outcome_data$outcome=outcomeName
			dat=harmonise_data(exposure_dat, outcome_data)
			dat=dat[dat$pval.outcome>5e-08,]
			#得到用于孟德尔随机化的工具变量
			snpTab=dat[dat$mr_keep=="TRUE",]
				
			if(nrow(snpTab)>=3){
				#孟德尔随机化分析
				mrResult=mr(dat)
				mrTab=generate_odds_ratios(mrResult)
				outTab=rbind(outTab, mrTab)
				pleioTab=mr_pleiotropy_test(dat)
				
				##########筛选阳性的基因,并对阳性的基因进行可视化#########
				#筛选IVW方法pvalue小于0.05的结果
				if((nrow(mrTab)>=3) && (mrResult$pval[3]<0.05)){
					#筛选五种方法OR方向一致的结果
					if(sum(mrTab$or>1)==nrow(mrTab) | sum(mrTab$or<1)==nrow(mrTab)){
						#多效性pvalue大于0.05
						if(as.numeric(pleioTab$pval)>0.05){
							#孟德尔随机化结果方向和差异分析方向一致
							if((geneRT[expoName,"logFC"] * mrResult$b[3])>0){
								#输出孟德尔随机化分析的结果
								write.csv(snpTab, file=paste0(i, ".table.SNP.csv"), row.names=F)
								write.csv(mrTab, file=paste0(i, ".table.MRresult.csv"), row.names=F)
								
								#MR-PRESSO异常值检测(偏倚的SNP)
								if(nrow(snpTab)>3){
									presso=run_mr_presso(dat)
									write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file=paste0(i, ".table.MR-PRESSO_Global.csv"))
									write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file=paste0(i, ".table.MR-PRESSO_Outlier.csv"))
								}
								
								#输出异质性检验的结果
								heterTab=mr_heterogeneity(dat)
								write.csv(heterTab, file=paste0(i, ".table.heterogeneity.csv"), row.names=F)
								
								#输出多效性检验的结果
								write.csv(pleioTab, file=paste0(i, ".table.pleiotropy.csv"), row.names=F)
								
								#绘制散点图
								pdf(file=paste0(i, ".scatter_plot.pdf"), width=7, height=6.5)
								p1=mr_scatter_plot(mrResult, dat)
								print(p1)
								dev.off()
								
								#森林图
								res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
								pdf(file=paste0(i, ".forest.pdf"), width=6.5, height=5)
								p2=mr_forest_plot(res_single)
								print(p2)
								dev.off()
								
								#漏斗图
								pdf(file=paste0(i, ".funnel_plot.pdf"), width=6.5, height=6)
								p3=mr_funnel_plot(singlesnp_results = res_single)
								print(p3)
								dev.off()
								
								#留一法敏感性分析
								pdf(file=paste0(i, ".leaveoneout.pdf"), width=6.5, height=5)
								p4=mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
								print(p4)
								dev.off()
							}
						}
					}
				}
			}
		}
	}
	#输出每个结局数据所有基因孟德尔随机化分析的结果
	write.csv(outTab, file=paste0(outcomeID, "__all.MRresult.csv"), row.names=F)
}

