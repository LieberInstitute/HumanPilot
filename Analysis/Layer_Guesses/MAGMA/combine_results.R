##

f = list.files("Results", pattern = "gsa.out",full=TRUE)
f = f[-5] # drop pgc3

names(f) = c("SZCD", "MDD", "ASD", "BPD")

gsaList = lapply(f, read.table, header=TRUE, as.is=TRUE,row.names=1,comment="#")

magmaTab = sapply(gsaList, "[[", "P")
rownames(magmaTab) = rownames(gsaList[[1]])

magmaTab = magmaTab[!grepl("Neuropil", rownames(magmaTab)),]
write.csv(magmaTab, file = "../suppTab_MAGMA_enrichment.csv")

magmaTabFdr = matrix(p.adjust(magmaTab, "fdr"), 
	nr = nrow(magmaTab),nc = ncol(magmaTab),
	dimnames = dimnames(magmaTab))

magmaTabBonf = matrix(p.adjust(magmaTab, "bonf"), 
	nr = nrow(magmaTab),nc = ncol(magmaTab),
	dimnames = dimnames(magmaTab))
