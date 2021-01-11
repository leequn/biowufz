# 设置环境变量
{
  rm(list = ls())
  setwd("/Users/liqun/tmp/")
}

# 加载数据以及前期处理
{
  # 读取数据
  tmpExpData <- read.table('tmp.txt',header = T)
  # 修改列名
  colnames(tmpExpData) <- c('Ensembl','symbol','oocyte_10days','oocyte_14days','GV_8weeks',
                         'MII_1','MII_2','PN5_1','PN5_2','twocell_early_1','twocell_early_2',
                         'twocell_late_1','twocell_late_2','fourcell_1','fourcell_2','eightcell_1',
                         'eightcell_2','ICM_1','ICM_2')
  # 创建表达矩阵，去重重复和空值
  tmpExpData <- tmpExpData[!duplicated(tmpExpData$symbol),]
  tmpExpData <- na.omit(tmpExpData)
  expData <- tmpExpData[,3:19]
  rownames(expData) <- tmpExpData$symbol
  # 数据处理和过滤，合子激活基因的筛选
  expData <- expData[expData$twocell_late_1>5,]
  expData <- expData[expData$MII_1<5,]
  expData <- expData[expData$twocell_late_2>5,]
  expData <- expData[expData$MII_2<5,]
}

# ZGA相关基因统计
{
  # 与ZGA相关基因
  ZGAgene <- expData
  # one-cell minor ZGA genes 
  onecellMinorZgagene <- expData[expData$PN5_1>2,]
  onecellMinorZgagene <- onecellMinorZgagene[onecellMinorZgagene$MII_1<0.5,]
  onecellMinorZgagene <- onecellMinorZgagene[onecellMinorZgagene$PN5_2>2,]
  onecellMinorZgagene <- onecellMinorZgagene[onecellMinorZgagene$MII_2<0.5,]
  # early two-cell minor ZGA genes
  twocellearlyMinorZgagene <- expData[expData$twocell_early_1>2,]
  twocellearlyMinorZgagene <- twocellearlyMinorZgagene[twocellearlyMinorZgagene$MII_1<0.5,]
  twocellearlyMinorZgagene <- twocellearlyMinorZgagene[twocellearlyMinorZgagene$twocell_early_2>2,]
  twocellearlyMinorZgagene <- twocellearlyMinorZgagene[twocellearlyMinorZgagene$MII_2<0.5,]
  
  num_zga=nrow(ZGAgene)
  num_onecell=nrow(onecellMinorZgagene)
  num_twocell=nrow(twocellearlyMinorZgagene)
  numbers = c(num_zga,num_onecell,num_twocell)
  
  barplot(numbers,
          main="Types of ZGA",
          col=c("#ED1C24","#22B14C","#FFC90E"),
          names.arg=c("All ZGA","one-cell minor ZGA","early two-cell minor ZGA")
  )
}

# 加载绘图包
{
  library(pheatmap)
  pheatmap(log(expData+0.0001),scale = "row",cluster_cols = FALSE,show_rownames=FALSE,main = "Heatmap of Mouse ZGA genes")
 
}

# 基因筛选
{
  # 根据差异表达结果（2cell晚期 vs MII期），筛选感兴趣的基因
  tmpData1 <- data.frame(expData$MII_1+expData$MII_2,expData$twocell_late_1+expData$twocell_late_2)
  rownames(tmpData1) <- rownames(expData)
  # 修改列名
  colnames(tmpData1) <- c("MII","twocell_late")
  tmpData1 <- tmpData1[tmpData1$twocell_late>2*tmpData1$MII,]
  tmpData1$genes <- rownames(tmpData1)
  
  # 筛选出在MII期、PN5期和2细胞早期表达变化不大的基因（FC<1.5）
  tmpData2 <- data.frame(expData$MII_1+expData$MII_2,expData$PN5_1+expData$PN5_2)
  rownames(tmpData2) <- rownames(expData)
  # 修改列名
  colnames(tmpData2) <- c("MII","PN5")
  tmpData2 <- tmpData2[tmpData2$PN5<1.5*tmpData2$MII,]
  tmpData2 <- tmpData2[tmpData2$PN5>2/3*tmpData2$MII,]
  tmpData2$genes <- rownames(tmpData2)
  
  tmpData3 <- data.frame(expData$twocell_early_1+expData$twocell_early_2,expData$PN5_1+expData$PN5_2)
  rownames(tmpData3) <- rownames(expData)
  # 修改列名
  colnames(tmpData3) <- c("twocell_early","PN5")
  tmpData3 <- tmpData3[tmpData3$PN5<1.5*tmpData3$twocell_early,]
  tmpData3 <- tmpData3[tmpData3$PN5>2/3*tmpData3$twocell_early,]
  tmpData3$genes <- rownames(tmpData3)
  
  # 合并文件
  tmpData_23 <-  merge(tmpData2,tmpData3,by='genes')
  finalData <-  merge(tmpData_23,tmpData1,by='genes')
  
  # 筛选目的基因，绘制热图
  tmpData4 <- expData
  tmpData4$genes <- rownames(expData)
  finalData1 <- merge(finalData,tmpData4,by='genes')
  # 提取相关信息
  finalData2 <- data.frame(finalData1$MII_1,finalData1$MII_2,finalData1$PN5_1,finalData1$PN5_2,finalData1$twocell_early_1,
                           finalData1$twocell_early_2,finalData1$twocell_late_1,finalData1$twocell_late_2)
  rownames(finalData2) <- finalData1$genes
  colnames(finalData2) <- c('MII_1',"MII_2","PN5_1","PN5_2","twocell_early_1","twocell_early_2","twocell_late_1","twocell_late_2")
  pheatmap(log2(finalData2+0.0001),scale = "row",cluster_cols = FALSE,show_rownames=FALSE,main = "Target ZGA genes")
}

# go以及kegg注释
{
  # 加载包
  library(org.Mm.eg.db)
  library(clusterProfiler)
  # 加载基因
  gene <- bitr(rownames(finalData2), fromType ="SYMBOL",
               toType =  "ENTREZID",
               OrgDb = org.Mm.eg.db)
  geneList <- bitr(tmpExpData$symbol, fromType ="SYMBOL",
                   toType =  "ENTREZID",
                   OrgDb = org.Mm.eg.db)
  # GO
  ego <- enrichGO(gene = gene$ENTREZID,
                  universe      = names(geneList$ENTREZID),
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  # KEGG
  kk <- enrichKEGG(gene         = gene$ENTREZID,
                   organism     = 'mmu',
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05)
  # 保存数据
  go = as.data.frame(ego)
  kegg = as.data.frame(kk)
  write.csv(go,'go.csv')
  write.csv(kegg,'kegg.csv')
}

# 绘图（GO和KEGG）
# ggplot2 setting
{
  ggplot_settings=function(legend=F){
    pp=theme(plot.title = element_text(hjust = 0.5,size=rel(1),face="bold"),
             panel.border = element_blank(), 
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), 
             panel.background = element_blank(),
             
             strip.text = element_text(size=rel(1),face="bold"),
             strip.background = element_blank(),
             
             axis.text.x =element_text(size=rel(0.8),angle=0),
             axis.text.y =element_text(size=rel(0.8),angle=0),
             axis.title=element_text(size=rel(0.8),face="bold"),
             axis.line = element_line(colour = "black", size = rel(0.7), linetype = "solid")
    )
    
    if(legend==F){
      pp=theme(legend.position = "none")+pp
    }
    return(pp)
  }
}
#GO
{
  # 加载包
  library(ggplot2)
  # 读取数据
  egoData <- read.csv('go.csv',header = T)
  keggData <- read.csv('kegg.csv',header = T)
  # 绘图
  ggplot(data = egoData)+
    geom_bar(aes(y=reorder(Description,Count),x=Count,fill=-log(pvalue)),stat='identity')+
    # Y轴为Term，以Count数排列，X轴为Count数；绘图函数里的stat参数表示对样本点做统计的方式，默认为identity，表示一个x对应一个y  
    scale_fill_gradient(expression(-log["10"](P.value)),low="grey",high="red")+
    # 设置图例
    ylab("")+
    # 设置背景
    xlab("Gene count")+ 
    # 设置背景
    ggplot_settings(T)+
    ggtitle("GO term of ZGA genes")
}
# kegg
{
  # 加载包
  #library(ggplot2)
  # 读取数据
  #egoData <- read.csv('go.csv',header = T)
  #keggData <- read.csv('kegg.csv',header = T)
  # 绘图
  ggplot(data = keggData)+
    geom_bar(aes(y=reorder(Description,Count),x=Count,fill=-log(pvalue)),stat='identity')+
    # Y轴为Term，以Count数排列，X轴为Count数；绘图函数里的stat参数表示对样本点做统计的方式，默认为identity，表示一个x对应一个y  
    scale_fill_gradient(expression(-log["10"](P.value)),low="grey",high="red")+
    # 设置图例
    ylab("")+
    # 设置背景
    xlab("Gene count")+ 
    # 设置背景
    ggplot_settings(T)+
    ggtitle("KEGG term of ZGA genes")
}



