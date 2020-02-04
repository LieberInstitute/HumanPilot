
## ----Libraries ------------------
library(SummarizedExperiment)
library(Matrix)


## load rse list
load("Human_DLPFC_Visium_processedData_rseList.rda")

## collect sample names and clusters
dfList = lapply(rseList, function(rse) {
	d = as.data.frame(colData(rse)[,c("sample_name", "Cluster")])
	d = d[!duplicated(d),]
	d = d[order(d$Cluster),]
	d
})
dd = do.call("rbind", dfList)
rownames(dd) = NULL
dd$Layer = ""

write.csv(dd, file = "guess_the_layer.csv",row.names=FALSE,quote=FALSE)
