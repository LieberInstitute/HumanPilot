## module load conda_R/3.6.x

## ----Libraries ------------------
library(parallel)
library(SummarizedExperiment)
library(Matrix)
library(RColorBrewer)
library(pdist) # for dist

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

## get distance w/in a sample
distList = mclapply(rseList, function(rse) {
	cat(".")
	u = assays(rse)$umis
	u = u[rowSums(u) > 0,]
	tu = t(u)
	dist(tu)
},mc.cores=4)
save(distList, file = "spotLevel_distances_withinImage.rda")

## get distance across reps
samplePairs = as.data.frame(combn(12,2))
distListPairs = mclapply(samplePairs, function(ii) {
	cat(".")
	u1 = assays(rseList[[ii[1]]])$umis
	u2 = assays(rseList[[ii[2]]])$umis
	tu1 = t(u1)
	tu2 = t(u2)
	as.matrix(pdist(tu1, tu2))
},mc.cores=4)
save(distListPairs, file = "spotLevel_distances_allPairs.rda")

