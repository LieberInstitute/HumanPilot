###
###
## module load conda_R/3.6.x

library(tidyverse)
library(ggplot2)
library(Matrix)
library(Rmisc)
library(ggforce)
library(cowplot)
library(RColorBrewer)
library(grid)
library(SummarizedExperiment)
library(jaffelab)
library(parallel)

## load rse list
load("Human_DLPFC_Visium_processedData_rseList.rda")

## add kmeans
kmFiles= list.files("/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/10X",
	pattern = "clustering", recur=TRUE, full=TRUE)
kmFiles = kmFiles[grep("kmeans", kmFiles)]

## get sample info
clusterDf = data.frame(SampleID = ss(kmFiles, "/", 8),
	K = as.numeric(ss(kmFiles, "_", 6)), 
	File = kmFiles, stringsAsFactors=FALSE)
clusterList = split(clusterDf, clusterDf$SampleID)
clusterList = clusterList[names(rseList)]

clList = lapply(clusterList, function(cc) {
	cc = cc[order(cc$K),]
	cl = lapply(cc$File, read.csv, 
		as.is=TRUE, row.names=1)
	clusters = do.call("cbind", cl)
	names(clusters) = paste0("k", cc$K)
	clusters
})

## join
rseList = mapply(function(rse, cl) {
	colData(rse) = cbind(colData(rse), cl)
	rse
}, rseList, clList)


## filter to expressed genes, lets be liberal
exprsMat = sapply(rseList, function(rse) {
	rowSums(assays(rse)$umis)
})

# hist(log2(rowSums(exprsMat)+1))
table(rowSums(exprsMat) > 0)
# liberal filter
exprsIndex = which(rowSums(exprsMat) > 0)
rseList = lapply(rseList, function(rse) rse[exprsIndex,])

## change rownames to symbols
rseList = lapply(rseList, function(rse) {
	rownames(rse) = rowData(rse)$gene_name
	rse
})
	
## counts per k
kTabList = lapply(2:10, function(k) {
	t(sapply(rseList, function(rse) {
		table(colData(rse)[,paste0("k", k)])
	}))
})


#####################
## start with more k, and then collapse to 7 
##   which assumes 6 layers and WM
#########

# sum up umis in each cluster
k=10
umiListComb = mclapply(rseList, function(rse) {
	cl = colData(rse)[,paste0("k", k)]
	u = assays(rse)$umis
	## prep counts
	lu = log2(u+1)
	su = scale(lu)

	## do kmeans on z-scale
	km = kmeans(t(su), centers = 7)
	
	## collapse once
	cIndexes = splitit(km$cluster)
	uSum = sapply(cIndexes, function(ii) rowSums(t(t(u[,ii]))))
	return(uSum)
	# collapse again
	# dd = dist(t(uSum))
	# hc = hclust(dd)
	# always for 7 groups
	# cl2 = cutree(hc,k=7)
	# cIndexes2 = splitit(cl2)
	# uSum2 = sapply(cIndexes2, function(ii) rowSums(t(t(uSum[,ii]))))
	# return(uSum2)
},mc.cores=6)

## corr within an image
lapply(umiListComb, function(x) cor(log2(x+1)))

u1 = umiListComb[[1]]
u2 = umiListComb[[2]]
cor(u1,u2)
u =  log2(umiListComb[[1]]+1)
cor(u)
plot(u[,1:2])
u["MOBP",]
### how do you line up groups?