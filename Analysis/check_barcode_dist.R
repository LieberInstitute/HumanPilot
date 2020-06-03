###
library(jaffelab)
library(Biostrings)

x = read.csv("../10X/151675/tissue_positions_list.txt",
	as.is=TRUE, header=FALSE)
	
x$barcode = ss(x$V1, "-")
bcs = DNAStringSet(x$barcode)

## hamming
dd_hamming = stringDist(bcs, method = "hamming")

dd_mat_hamming = as.matrix(dd_hamming)
dd_mat_hamming[upper.tri(dd_mat_hamming, diag=TRUE)] = NA

min_dist_hamming = apply(dd_mat_hamming, 1, min, na.rm=TRUE)
table(min_dist_hamming)

## 
dd = stringDist(bcs, method = "levenshtein")

dd_mat = as.matrix(dd)
dd_mat[upper.tri(dd_mat, diag=TRUE)] = NA

min_dist = apply(dd_mat, 1, min, na.rm=TRUE)
table(min_dist)

checks = which(dd_mat == 3, arr.ind=TRUE)
checks = as.data.frame(checks)
checks$bc1 = x$barcode[checks$row]
checks$bc2 = x$barcode[checks$col]
rownames(checks)= NULL
colnames(checks)[1:2] = paste0("row", 1:2)

checks