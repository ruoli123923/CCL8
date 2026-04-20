#设置工作目录
setwd("~/")

#读入,因后续有标准化操作，故不需标准化
GSE217314.exp=read.table("predata/GSE217314/GSE217314.txt", header=T, sep="\t", check.names=F,row.names = 1)
GSE224570.exp=read.table("predata/GSE224570/GSE224570.txt", header=T, sep="\t", check.names=F,row.names = 1)
GSE206932.exp=read.table("predata/GSE206932/GSE206932.txt", header=T, sep="\t", check.names=F,row.names = 1)

gene_set_200 <- read.table(
  "Ccl8_DEG_top200.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)
gene_set_200_vec <- gene_set_200$V1

#获取共有基因
samegene = Reduce(intersect,list(rownames(GSE206932.exp),rownames(GSE224570.exp),rownames(GSE217314.exp)))
samegene <- intersect(gene_set_200_vec, samegene)

#分别提取
GSE206932.exp = GSE206932.exp[samegene,]
GSE224570.exp = GSE224570.exp[samegene,]
GSE217314.exp = GSE217314.exp[samegene,]

#合并，导出
Test <- cbind(GSE206932.exp, GSE224570.exp)
Train_expr <- GSE217314.exp

#加载包
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
library(limma)
library(ggpubr)

#加载模型训练以及模型评估的脚本
source("ML.R")
Train_expr <- GSE217314.exp

common_genes <- Reduce(intersect, list(
  rownames(GSE206932.exp),
  rownames(GSE224570.exp)
))

Test_expr <- cbind(GSE206932.exp, GSE224570.exp)

#读入训练集表达矩阵，及分组信息（结局事件）
Train_class <- read.table(
  "predata/Training_class.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(Train_class) <- Train_class$sample

#读入训练集表达矩阵，及分组信息（结局事件）
Test_class <- read.table("predata/Testing_class.txt", header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)

#提取共有样本
comsam <- intersect(rownames(Test_class), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]
Test_class <- Test_class[comsam,,drop = F]

#提取相同基因
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[common_genes,])
Test_expr <- t(Test_expr[common_genes,])

#标准化
Train_set <- t(scale(t(Train_expr), center = TRUE, scale = TRUE))
Test_set  <- t(scale(t(Test_expr), center = TRUE, scale = TRUE))
#去除NA然后重新保留样本
Train_set <- Train_set[complete.cases(Train_set), ]
Train_class <- Train_class[complete.cases(Train_set), ]
common_samples <- intersect(rownames(Train_set), Train_class$sample)
Train_set <- Train_set[common_samples, ]
Train_class <- Train_class[Train_class$sample %in% common_samples, ]

#获取模型名的列表
methods <- read.xlsx("methods.xlsx", startRow = 2)
methods <- methods$Model
methods <- gsub("-| ", "", methods)

#设置，结局事件的列名，最小的纳入基因数
classVar = "outcome"
min.selected.var = 3 

#预训练，先运行变量筛选的算法
preTrain.var <- list()
set.seed(seed = 777)

Variable = colnames(Train_set)
preTrain.method =  strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1]) 
preTrain.method = unique(unlist(preTrain.method)) 
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, 
                                 Train_set = Train_set, 
                                 Train_label = Train_class, 
                                 #Variable(筛选变量)和Model(获取模型)
                                 mode = "Variable", 
                                 classVar = classVar)
}

preTrain.var[["simple"]] <- colnames(Train_set)

#模型训练
model <- list()
set.seed(seed = 777)
Train_set_bk = Train_set

for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  
  method_name = method
  method <- strsplit(method, "\\+")[[1]]
  if (length(method) == 1) method <- c("simple", method) 
  
  Variable = preTrain.var[[method[1]]]
  Train_set = Train_set_bk[, Variable]
  Train_label = data.frame(outcome = Train_class$outcome)
  
  model[[method_name]] <- RunML(
    method = method[2],
    Train_set = Train_set,
    Train_label = Train_label,
    mode = "Model",
    classVar = "outcome"
  )
  
  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}

Train_set = Train_set_bk

#保存
saveRDS(model, "model.rds") 

#读取
model <- readRDS( "model.rds")
methodsValid <- names(model)

#计算样本风险评分
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set))
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(data.frame(ID=rownames(RS_mat),RS_mat), "RS_mat.txt",sep = "\t", row.names = F, quote = F)

#预测分类
Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(fit = model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set))
}

Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))

write.table(data.frame(ID=rownames(Class_mat),Class_mat),"Class_mat.txt", sep = "\t", row.names = F, quote = F)

#获取模型中的基因列表
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

write.table(fea_df, "fea_df.txt",sep = "\t", row.names = F, col.names = T, quote = F)

#计算各模型的C-index
Test_set <- as.data.frame(Test_set)
Train_set <- as.data.frame(Train_set)
AUC_list <- list()
for (method in methodsValid){
  AUC_list[[method]] <- RunEval(fit = model[[method]],
                                Test_set = Test_set,
                                Test_label = Test_class,
                                #注意修改
                                Train_name = "Training", 
                                #计算A平均AUC值，是否纳入训练集
                                Train_set = Train_set,
                                Train_label = Train_class,
                                # Train_set = NULL,
                                # Train_label = NULL,
                                cohortVar = "cohort",
                                classVar = classVar)
}
AUC_mat <- do.call(rbind, AUC_list)
write.table(data.frame(ID=rownames(AUC_mat),AUC_mat), "AUC_mat.txt",sep = "\t", row.names = F, col.names = T, quote = F)

#热图
AUC_mat <- read.table("AUC_mat.txt",sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_AUC <- apply(AUC_mat, 1, mean) 
avg_AUC <- sort(avg_AUC, decreasing = T) 
AUC_mat <- AUC_mat[names(avg_AUC), ] 
fea_sel <- fea_list[[rownames(AUC_mat)[1]]]
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3)) 
#颜色
if(ncol(AUC_mat) < 3) { 
  CohortCol <- c("red","blue")
} else { 
  CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired") 
}
names(CohortCol) <- colnames(AUC_mat)
cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(AUC_mat, 
                    avg_AUC, 
                    CohortCol, "steelblue", 
                    cellwidth = cellwidth, cellheight = cellheight, 
                    cluster_columns = F, cluster_rows = F) 
pdf("AUC.pdf", width = cellwidth * ncol(AUC_mat) + 7, height = cellheight * nrow(AUC_mat) * 0.45)
draw(hm)
invisible(dev.off())

####diff####

#首先需读入499-548行
#模型中基因的差异分析
#读入模型中的基因
model.gene=read.table("fea_df.txt", header=T, sep="\t", check.names=F)
#最优算法的名字
best.model = "Ridge"
best.model.gene = subset(model.gene,algorithm == best.model)[,1]

#训练集
#表达矩阵
Train = t(Train_expr[,best.model.gene])

#注意修改数据集名字
genesets = "Traning"

#分组信息
con = subset(Train_class,outcome == 0)
treat = subset(Train_class,outcome == 1)
conData=Train[,rownames(con)]
treatData=Train[,rownames(treat)]
Train=cbind(conData, treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)
#修改分组信息
Type=c(rep("Young",conNum), rep("Aged",treatNum))
my_comparisons=list()
my_comparisons[[1]]=levels(factor(Type))

#差异分析
newGeneLists=c()
outTab=data.frame()

for(i in row.names(Train)){
  rt1=data.frame(expression=Train[i,], Type=Type)
  
  #对差异基因进行可视化，绘制箱线图
  boxplot=ggboxplot(rt1, x="Type", y="expression", color="Type",
                    xlab="",
                    ylab=paste(i, "expression"),
                    legend.title="",
                    palette = c("#00AF50","#F5B700"),
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons,method = "t.test")
  pdf(file=paste0("diff/",genesets,".diff.",i,".pdf"), width=5, height=4.5)
  print(boxplot)
  dev.off()
}





#测试集
#表达矩阵
Test = t(Test_expr[,best.model.gene])

#注意修改数据集名字
unique(Test_class$cohort)
genesets = "GSE206932"

#分组信息
con = subset(Test_class,outcome == 0 &  cohort == genesets)[,2,drop=F]
treat = subset(Test_class,outcome == 1&  cohort == genesets)[,2,drop=F]
conData=Test[,rownames(con)]
treatData=Test[,rownames(treat)]
Test=cbind(conData, treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

#修改分组信息
Type=c(rep("Young",conNum), rep("Aged",treatNum))
my_comparisons=list()
my_comparisons[[1]]=levels(factor(Type))

#差异分析
newGeneLists=c()
outTab=data.frame()

for(i in row.names(Test)){
  rt1=data.frame(expression=Test[i,], Type=Type)
  
  #对差异基因进行可视化，绘制箱线图
  boxplot=ggboxplot(rt1, x="Type", y="expression", color="Type",
                    xlab="",
                    ylab=paste(i, "expression"),
                    legend.title="",
                    palette = c("#00AF50","#F5B700"),
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons,method = "t.test")
  pdf(file=paste0("diff/",genesets,".diff.",i,".pdf"), width=5, height=4.5)
  print(boxplot)
  dev.off()
}


####ROC####

#首先需读入499-548行
#模型中基因的ROC曲线
#读入模型中的基因
model.gene=read.table("fea_df.txt20417", header=T, sep="\t", check.names=F)
#最优算法的名字
best.model = "Ridge"
best.model.gene = subset(model.gene,algorithm == best.model)[,1]

#训练集
#表达矩阵
Train = t(Train_set[,best.model.gene])

#注意修改数据集名字
genesets = "Traning"

#分组信息
con = subset(Train_class,outcome == 0)
treat = subset(Train_class,outcome == 1)
conData=Train[,rownames(con)]
treatData=Train[,rownames(treat)]
Train=cbind(conData, treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)
y=c(rep(0,conNum), rep(1,treatNum))

#对交集基因进行循环，绘制ROC曲线
for(x in best.model.gene){
  #绘制ROC曲线
  roc1=roc(y, as.numeric(Train[x,]))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("ROC/",genesets,".ROC.",x,".pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=x)
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}

#测试集
#表达矩阵
Test = t(Test_set[,best.model.gene])

#注意修改数据集名字
unique(Test_class$cohort)
genesets = "GSE224570"

#分组信息
con = subset(Test_class,outcome == 0 &  cohort == genesets)[,2,drop=F]
treat = subset(Test_class,outcome == 1&  cohort == genesets)[,2,drop=F]
conData=Test[,rownames(con)]
treatData=Test[,rownames(treat)]
Test=cbind(conData, treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)
y=c(rep(0,conNum), rep(1,treatNum))

#对交集基因进行循环，绘制ROC曲线
for(x in best.model.gene){
  #绘制ROC曲线
  roc1=roc(y, as.numeric(Test[x,]))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  pdf(file=paste0("ROC/",genesets,".ROC.",x,".pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=x)
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}


