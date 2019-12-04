## module load conda_R/3.6.x

## ----Libraries ------------------
library(parallel)
library(SummarizedExperiment)
library(Matrix)
library(RColorBrewer)
library(jaffelab) 
library(edgeR)

## load rse list
load("Human_DLPFC_Visium_processedData_rseList.rda")

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
	
#### k means
Ks = 6:12
names(Ks) = Ks
kList = lapply(Ks, function(k) {
	cat(k)
	mclapply(rseList, function(rse) {
		cat(".")
		## prep counts
		u = assays(rse)$umis
		lu = log2(u+1)
		su = scale(lu)

		## do kmeans on z-scale
		km = kmeans(t(su), centers = k)
		
		## get DE by cluster
		cIndexes = splitit(km$cluster)
		tMat = sapply(cIndexes, function(ii) {
			x = rep(0,ncol(rse))
			x[ii] = 1
			f = lmFit(su, model.matrix(~x))
			f$coef[,2]
		})
		
		list(kmeans = km, tMat = tMat)
	},mc.cores=4)
})

save(kList, file = "kmeans_k6to12_listOutput_withTstat.rda")