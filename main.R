#filtering parameters
minCellNum           = 3                                                # filtering, remove genes expressed in fewer than minCellNum cells
expressed_cutoff     = 1                                                # filtering, for raw counts

#kmeans parameters
K.max                = 7                                               # if using the gap statistic, highest k that should be considered

#Gini parameters
gini.bi              = 0                                                # fitting, default is 0, for qPCR data, set as 1. log2.expr.cutoffl    = 0                                                # cutoff for range of gene expression   
log2.expr.cutoffh    = 20                                               # cutoff for range of gene expression 
log2.expr.cutoffl    = 0                                                # cutoff for range of gene expression
span                 = 0.9                                              # parameter for LOESS fitting
outlier_remove       = 0.75                                             # parameter for LOESS fitting
Gini.pvalue_cutoff   = 0.0001                                           # fitting, Pvalue, control how many Gini genes chosen
Gamma                = 0.99                                             # parameter for clustering

#ProgClust parameters

crit=3.5                                                               # consistent with ros in the paper, criteria for identifying rare cells.. 
o_crit=0.2                                                             # consistent with os in the paper, criteria for identifying rare cells.. 
k_crit =0.12                                                           # control the number of clusters. 
filt_c=1                                                               # filter cells with low expression levels.

# packages
library(speccalt)
library(apcluster)

#load packages and functions
source("F:/LH/code/GiniClust2_packages.R")
source("F:/LH/code/ISClust_GiniClust2_functions.R")

#Simulation data

ExprM.RawCounts <- read.table("F:/LH/GiniClust2/Proj/Simulation/data/Data_2000_1000_10_6_4_3.xls", sep="\t", head=TRUE, row.names=1)

# # mouse embryonic stem cells
# ExprM.RawCounts <- read.delim("F:/LH/GiniClust2/Proj/breast2/GSM1599498_ES_d4_LIFminus.csv", sep=",", head=F)
# ExprM.RawCounts.b<-ExprM.RawCounts
# ExprM.RawCounts<-ExprM.RawCounts2
# #raw data
# title=c("Symbol");
# for(i in 2:ncol(ExprM.RawCounts)){
#   title=c(title,paste(exprimentID, ".Cell_",i-1,sep=""))
# };
# colnames(ExprM.RawCounts)=title
# rownames(ExprM.RawCounts)=ExprM.RawCounts[,1]
# ExprM.RawCounts= ExprM.RawCounts[,-1]


dim(ExprM.RawCounts)

#filtering
ExpressedinCell_per_gene=apply(ExprM.RawCounts,1,function(x) length(x[x > expressed_cutoff ]))
ExprM.RawCounts.filter = ExprM.RawCounts[ExpressedinCell_per_gene >= minCellNum,]

#clustering
finalCluster<-igclust(ExprM.RawCounts.filter,deep=1)

label<-c(rep(1,2000),rep(2,1000),rep(3,10),rep(4,6),rep(5,4),rep(6,3))
NMI<-wnmi(finalCluster,label)
purity<-wpurity(finalCluster,label)