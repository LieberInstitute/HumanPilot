###

pd = read.csv("he_SraRunTable.txt",as.is=TRUE)

## number of unique variables?
uniqueIndex = which(sapply(pd, function(x) length(unique(x))) > 1)
apply(pd[,uniqueIndex], 2, table)
table(pd$Sample_Name, pd$tissue)

fqPath = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/he_layers/FASTQ/"

man = data.frame(leftRead = paste0(fqPath, pd$Run, "_1.fastq.gz"),
	leftMmd5 = 0, sampleID = pd$Run, stringsAsFactors=FALSE)
all(file.exists(man$leftRead))
all(file.exists(man$rightRead))

write.table(man, file="samples.manifest", sep="\t",
	quote=FALSE, col.names=FALSE, row.names=FALSE)