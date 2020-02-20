###

pd = read.delim("hafner_SraRunTable.txt",as.is=TRUE)
table(pd$Sample_Name, pd$tissue)

fqPath = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/hafner_vglut/FASTQ/"

man = data.frame(leftRead = paste0(fqPath, pd$Run, "_1.fastq.gz"),
	leftMmd5 = 0, sampleID = pd$Run, stringsAsFactors=FALSE)
all(file.exists(man$leftRead))
all(file.exists(man$rightRead))

write.table(man, file="samples.manifest", sep="\t",
	quote=FALSE, col.names=FALSE, row.names=FALSE)