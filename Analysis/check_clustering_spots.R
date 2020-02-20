## module load conda_R/3.6.x

## ----Libraries ------------------
library(parallel)
library(SummarizedExperiment)
library(Matrix)
library(RColorBrewer)
library(pdist) # for dist
library(
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

####################
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

############# read back in ##########
load("spotLevel_distances_withinImage.rda")
load("spotLevel_distances_allPairs.rda")

mean_dist_pairs = sapply(distListPairs, mean)
dat = data.frame(pairs = t(samplePairs), meanDist = mean_dist_pairs)
colnames(dat)[1:2] = paste0("Sample", 1:2)

dat$Sample1 = names(rseList)[dat$Sample1]
dat$Sample2 = names(rseList)[dat$Sample2]

## add more metrics
load('Human_DLPFC_Visium_processedData_sce_scran.Rdata')
metrics = colData(sce)
metrics = metrics[!duplicated(metrics$sample_name),]

dat$SubjPos1 = metrics$subject_position[match(dat$Sample1, metrics$sample_name)]
dat$SubjPos2 = metrics$subject_position[match(dat$Sample2, metrics$sample_name)]
dat$Subj1 = metrics$subject[match(dat$Sample1, metrics$sample_name)]
dat$Subj2 = metrics$subject[match(dat$Sample2, metrics$sample_name)]

dat$spatialDup = dat$SubjPos1 == dat$SubjPos2
dat$samePerson = dat$Subj1 == dat$Subj2
dat$lab = factor(paste0(dat$spatialDup,":", dat$samePerson))

dat[dat$lab == "TRUE:TRUE",]
boxplot(meanDist ~ lab, data=dat,outline=FALSE)
points(meanDist ~ jitter(as.numeric(lab),0.15), data=dat,pch=21,bg="grey")

summary(lm(meanDist ~ spatialDup + samePerson, data=dat))