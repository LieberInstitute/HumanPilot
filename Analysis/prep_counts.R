library('derfinder')
library('BiocParallel')
library('Biostrings')
library('GenomicRanges')
library('GenomicFeatures')
library('org.Hs.eg.db')
library('biomaRt')
library('jaffelab')
library('getopt')
library('rafalib')
library('devtools')
library('SummarizedExperiment')
library('plyr')
library('rtracklayer')

## read in pheno
manifest <- read.table('/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/sample_level_layer_map.tsv', 
	sep = '\t',header = FALSE, stringsAsFactors = FALSE)

colnames(manifest) = c("SlideID", "Layer", "bam")
manifest$SampleID = paste0(manifest$SlideID, "_", manifest$Layer)

## load ref
gtf = "/dcl01/ajaffe/data/lab/singleCell/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf"
gencodeGTF = import(gtf)
gencodeGENES = mcols(gencodeGTF)[which(gencodeGTF$type=="gene"),
	c("gene_id","type","gene_biotype","gene_name")]
rownames(gencodeGENES) = gencodeGENES$gene_id

gencodeEXONS = as.data.frame(gencodeGTF)[which(gencodeGTF$type=="exon"),
	c("seqnames","start","end","gene_id","exon_id")]
names(gencodeEXONS) = c("Chr","Start","End","gene_id","exon_gencodeID")

### exons in PAR regions
par_y = grep("PAR_Y",gencodeEXONS$gene_id)
gencodeEXONS$exon_gencodeID[par_y] = paste0(gencodeEXONS$exon_gencodeID[par_y],"_PAR_Y")
gencodeEXONS = gencodeEXONS[,-4]

###############
### gene counts
geneFn <- paste0("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/",
	manifest$SlideID, "/Layers/", manifest$SampleID, ".genes.counts")
names(geneFn) = manifest$SampleID
stopifnot(all(file.exists(geneFn)))

### read in annotation ##
geneMap = read.delim(geneFn[1], skip=1, as.is=TRUE)[,1:6]

## organize gene map
geneMap$Chr = ss(geneMap$Chr, ";")
geneMap$Start = as.numeric(ss(geneMap$Start, ";"))
tmp = strsplit(geneMap$End, ";")
geneMap$End = as.numeric(sapply(tmp, function(x) x[length(x)]))
geneMap$Strand = ss(geneMap$Strand, ";")
rownames(geneMap) = geneMap$Geneid
geneMap$gencodeID = geneMap$Geneid
geneMap$ensemblID = ss(geneMap$Geneid, "\\.")
geneMap$Geneid = NULL
geneMap$gene_type = gencodeGENES[geneMap$gencodeID,"gene_type"]
geneMap$Symbol = gencodeGENES[geneMap$gencodeID,"gene_name"]

######### biomart 
if (opt$organism=="hg19") {
	# VERSION 75, GRCh37.p13
	ensembl = useMart("ENSEMBL_MART_ENSEMBL", 
		dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
	sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), 
		values=geneMap$ensemblID, mart=ensembl)
} else if (opt$organism=="hg38") {
	# VERSION 85, GRCh38.p7
	ensembl = useMart("ENSEMBL_MART_ENSEMBL",  
		dataset="hsapiens_gene_ensembl", host="jul2016.archive.ensembl.org")
	sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), 
			values=geneMap$ensemblID, mart=ensembl)
}
#########

# geneMap$Symbol = sym$hgnc_symbol[match(geneMap$ensemblID, sym$ensembl_gene_id)]
geneMap$EntrezID = sym$entrezgene[match(geneMap$ensemblID, sym$ensembl_gene_id)]

## counts
geneCountList = mclapply(geneFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=opt$cores)
geneCounts = do.call("cbind", geneCountList)
rownames(geneCounts) = rownames(geneMap)
geneCounts = geneCounts[,metrics$SAMPLE_ID] # put in order

# number of reads assigned
geneStatList = lapply(paste0(geneFn, ".summary"), 
	read.delim,row.names=1)
geneStats = do.call("cbind", geneStatList)
colnames(geneStats) = metrics$SAMPLE_ID
metrics$totalAssignedGene = as.numeric(geneStats[1,] / colSums(geneStats))
## Add all the other stats from featureCounts at the gene level
geneStats_t <- t(geneStats)
colnames(geneStats_t) <- paste0('gene_', colnames(geneStats_t))
metrics <- cbind(metrics, geneStats_t)

# rna Rate
metrics$rRNA_rate = colSums(geneCounts[which(geneMap$gene_type == "rRNA"),])/colSums(geneCounts)


# make RPKM
bg = matrix(rep(as.numeric(geneStats["Assigned",])), nc = nrow(metrics), 
	nr = nrow(geneCounts),	byrow=TRUE)
widG = matrix(rep(geneMap$Length), nr = nrow(geneCounts), 
	nc = nrow(metrics),	byrow=FALSE)
geneRpkm = geneCounts/(widG/1000)/(bg/1e6)

## save metrics
write.csv(metrics, file = file.path(opt$maindir,
    paste0('read_and_alignment_metrics_', opt$experiment, '_', opt$prefix,
    '.csv')))


###############
### exon counts
exonFn <- file.path(opt$maindir, 'Counts', 'exon', paste0(metrics$SAMPLE_ID, filename, '_Exons.counts'))
names(exonFn) = metrics$SAMPLE_ID
stopifnot(all(file.exists(exonFn)))

### read in annotation ##
exonMap = read.delim(exonFn[1], skip=1, as.is=TRUE)[,1:6]
exonMap$gencodeID = exonMap$Geneid
exonMap$ensemblID = ss(exonMap$Geneid, "\\.")
rownames(exonMap) = paste0("e", rownames(exonMap))
exonMap$Geneid = NULL
exonMap$gene_type = gencodeGENES[exonMap$gencodeID,"gene_type"]
exonMap$Symbol = gencodeGENES[exonMap$gencodeID,"gene_name"]

# exonMap$Symbol = sym$hgnc_symbol[match(exonMap$ensemblID, sym$ensembl_gene_id)]
exonMap$EntrezID = sym$entrezgene[match(exonMap$ensemblID, sym$ensembl_gene_id)]

## add gencode exon id
exonMap = join(exonMap, gencodeEXONS, type="left", match="first")

## counts
exonCountList = mclapply(exonFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=opt$cores)
exonCounts = do.call("cbind", exonCountList)
rownames(exonCounts) = rownames(exonMap)
exonCounts = exonCounts[,metrics$SAMPLE_ID] # put in order

## remove duplicated
eMap = GRanges(exonMap$Chr, IRanges(exonMap$Start, exonMap$End))

## drop runthrough exons with duplicated exons
i = grepl('-', exonMap$Symbol)
j = countOverlaps(eMap[i], eMap[!i], type = 'equal') > 0
dropIndex = which(i)[j]
exonCounts = exonCounts[-dropIndex,]
exonMap = exonMap[-dropIndex,]
eMap = eMap[-dropIndex,]

## drop duplicated exons
keepIndex = which(!duplicated(eMap))
exonCounts = exonCounts[keepIndex,]
exonMap = exonMap[keepIndex,]

## change rownames
exonMap$exon_libdID = rownames(exonMap)
# rownames(exonMap) = rownames(exonCounts) = exonMap$exon_gencodeID # not always unique

# number of reads assigned
exonStatList = lapply(paste0(exonFn, ".summary"), 
                      read.delim,row.names=1)
exonStats = do.call("cbind", exonStatList)
colnames(exonStats) = metrics$SAMPLE_ID

## make RPKM
bgE = matrix(rep(colSums(exonCounts)), nc = nrow(metrics), 
	nr = nrow(exonCounts),	byrow=TRUE)
widE = matrix(rep(exonMap$Length), nr = nrow(exonCounts), 
	nc = nrow(metrics),	byrow=FALSE)
exonRpkm = exonCounts/(widE/1000)/(bgE/1e6)


############################
### add transcript maps ####
if (opt$organism == "hg19") { 
	load(file.path(RDIR, "feature_to_Tx_hg19_gencode_v25lift37.rda"))
} else if (opt$organism == "hg38") { 
	load(file.path(RDIR, "feature_to_Tx_hg38_gencode_v25.rda")) 
	load(file.path(RDIR, "exonMaps_by_coord_hg38_gencode_v25.rda"))
}

## gene annotation
geneMap$Class = "InGen"
geneMap$meanExprs = rowMeans(geneRpkm)
mmTx = match(geneMap$gencodeID, names(allTx))
tx = CharacterList(vector("list", nrow(geneMap)))
tx[!is.na(mmTx)] = allTx[mmTx[!is.na(mmTx)]]
geneMap$NumTx = elementNROWS(tx)
geneMap$gencodeTx = sapply(tx,paste0,collapse=";")

## exon annotation
exonMap$Class = "InGen"
exonMap$meanExprs = rowMeans(exonRpkm)
exonMap$coord = paste0(exonMap$Chr,":",exonMap$Start,"-",exonMap$End,"(",exonMap$Strand,")")
if (opt$organism == "hg19") {
	mmTx = match(exonMap$exon_libdID, names(allTx))
	tx = CharacterList(vector("list", nrow(exonMap)))
	tx[!is.na(mmTx)] = allTx[mmTx[!is.na(mmTx)]]
	exonMap$NumTx = elementNROWS(tx)
	exonMap$gencodeTx = sapply(tx,paste0,collapse=";")
	
} else if (opt$organism == "hg38") { 
	exonMap = exonMap[,-which(colnames(exonMap) %in% c("exon_gencodeID","exon_libdID"))]

	mmENSE = match(exonMap$coord, names(coordToENSE))
	ENSE = CharacterList(vector("list", nrow(exonMap)))
	ENSE[!is.na(mmENSE)] = coordToENSE[mmENSE[!is.na(mmENSE)]]
	exonMap$NumENSE = elementNROWS(ENSE)
	exonMap$exon_gencodeID = sapply(ENSE,paste0,collapse=";")
	
	mmLIBD = match(exonMap$coord, names(coordToEid))
	libdID = CharacterList(vector("list", nrow(exonMap)))
	libdID[!is.na(mmLIBD)] = coordToEid[mmLIBD[!is.na(mmLIBD)]]
	exonMap$NumLIBD = elementNROWS(libdID)
	exonMap$exon_libdID = sapply(libdID,paste0,collapse=";")
	
	mmTx = match(exonMap$coord, names(coordToTX))
	tx = CharacterList(vector("list", nrow(exonMap)))
	tx[!is.na(mmTx)] = coordToTX[mmTx[!is.na(mmTx)]]
	exonMap$NumTx = elementNROWS(tx)
	exonMap$gencodeTx = sapply(tx,paste0,collapse=";")
}


## Create gene,exon RangedSummarizedExperiment objects

gr_genes <- GRanges(seqnames = geneMap$Chr,
    IRanges(geneMap$Start, geneMap$End), strand = geneMap$Strand)
names(gr_genes) <- rownames(geneMap)
mcols(gr_genes) <- DataFrame(geneMap[, - which(colnames(geneMap) %in%
    c('Chr', 'Start', 'End', 'Strand'))])

rse_gene <- SummarizedExperiment(assays = list('counts' = geneCounts),
    rowRanges = gr_genes, colData = metrics)
save(rse_gene, file = paste0('rse_gene_', EXPNAME, '_n', N, '.Rdata'))

gr_exons <- GRanges(seqnames = exonMap$Chr,
    IRanges(exonMap$Start, exonMap$End), strand = exonMap$Strand)
names(gr_exons) <- rownames(exonMap)
mcols(gr_exons) <- DataFrame(exonMap[, - which(colnames(exonMap) %in%
    c('Chr', 'Start', 'End', 'Strand'))])
    
rse_exon <- SummarizedExperiment(assays = list('counts' = exonCounts),
    rowRanges = gr_exons, colData = metrics)
save(rse_exon, file = paste0('rse_exon_', EXPNAME, '_n', N, '.Rdata'))



###################
##### junctions

## import theJunctions annotation
if (opt$organism == "hg19") { 
	load(file.path(RDIR, "junction_annotation_hg19_gencode_v25lift37.rda"))
} else if (opt$organism == "hg38") { 
	load(file.path(RDIR, "junction_annotation_hg38_gencode_v25.rda"))
}

## via primary alignments only
junctionFiles <- file.path(opt$maindir, 'Counts', 'junction', paste0(metrics$SAMPLE_ID, '_junctions_primaryOnly_regtools.count'))
stopifnot(all(file.exists(junctionFiles))) #  TRUE

if (opt$stranded %in% c('forward', 'reverse')) {
	juncCounts = junctionCount(junctionFiles, metrics$SAMPLE_ID,
		output = "Count", maxCores=opt$cores,strandSpecific=TRUE)
} else {
	juncCounts = junctionCount(junctionFiles, metrics$SAMPLE_ID,
		output = "Count", maxCores=opt$cores,strandSpecific=FALSE)
}
## filter junction counts - drop jxns in <1% of samples
n = max(1, floor(N/100))
jCountsLogical = DataFrame(sapply(juncCounts$countDF, function(x) x > 0))
jIndex = which(rowSums(as.data.frame(jCountsLogical)) >= n) 
juncCounts = lapply(juncCounts, function(x) x[jIndex,])


############ anno/jMap
anno = juncCounts$anno
# seqlevels(anno, force=TRUE) = paste0("chr", c(1:22,"X","Y","M"))
seqlevels(anno, pruning.mode="coarse") = paste0("chr", c(1:22,"X","Y","M"))  # for updated R 3.5

## add additional annotation
anno$inGencode = countOverlaps(anno, theJunctions, type="equal") > 0
anno$inGencodeStart = countOverlaps(anno, theJunctions, type="start") > 0
anno$inGencodeEnd = countOverlaps(anno, theJunctions, type="end") > 0

oo = findOverlaps(anno, theJunctions, type="equal")
anno$gencodeGeneID = NA
anno$gencodeGeneID[queryHits(oo)] = as.character(theJunctions$gencodeID[subjectHits(oo)])
anno$ensemblID = ss(anno$gencodeGeneID, "\\.")
anno$Symbol = NA
anno$Symbol[queryHits(oo)] = theJunctions$symbol[subjectHits(oo)]
anno$gencodeStrand = NA
anno$gencodeStrand[queryHits(oo)] = as.character(strand(theJunctions)[subjectHits(oo)])
anno$gencodeTx = CharacterList(vector("list", length(anno)))
anno$gencodeTx[queryHits(oo)] = theJunctions$tx[subjectHits(oo)]
anno$numTx = elementNROWS(anno$gencodeTx)

## junction code
anno$code = ifelse(anno$inGencode, "InGen", 
	ifelse(anno$inGencodeStart & anno$inGencodeEnd, "ExonSkip",
	ifelse(anno$inGencodeStart | anno$inGencodeEnd, "AltStartEnd", "Novel")))

## b/w exons and junctions
exonGR = GRanges( exonMap$Chr,	IRanges(exonMap$Start, exonMap$End))
anno$startExon = match(paste0(seqnames(anno),":",start(anno)-1), 
	paste0(seqnames(exonGR), ":", end(exonGR)))
anno$endExon = match(paste0(seqnames(anno),":",end(anno)+1),
	paste0(seqnames(exonGR), ":", start(exonGR)))
g = data.frame(leftGene = exonMap$gencodeID[anno$startExon],
	rightGene = exonMap$gencodeID[anno$endExon],
	leftGeneSym = exonMap$Symbol[anno$startExon],
	rightGeneSym = exonMap$Symbol[anno$endExon],
	stringsAsFactors=FALSE)
g$newGene = NA
g$newGeneSym = NA
g$newGene[which(g$leftGene==g$rightGene)] = 
	g$leftGene[which(g$leftGene==g$rightGene)] 
g$newGeneSym[which(g$leftGene==g$rightGene)] = 
	g$leftGeneSym[which(g$leftGene==g$rightGene)] 
g$newGene[which(g$leftGene!=g$rightGene)] = 
	paste0(g$leftGene,"-",g$rightGene)[which(g$leftGene!=g$rightGene)] 
g$newGeneSym[which(g$leftGene!=g$rightGene)] = 
	paste0(g$leftGeneSym,"-",g$rightGeneSym)[which(g$leftGene!=g$rightGene)] 
g$newGene[which(is.na(g$newGene) & is.na(g$leftGene))] = 
	g$rightGene[which(is.na(g$newGene) & is.na(g$leftGene))] 
g$newGene[which(is.na(g$newGene) & is.na(g$rightGene))] = 
	g$leftGene[which(is.na(g$newGene) & is.na(g$rightGene))] 
g$newGeneSym[which(is.na(g$newGeneSym) & is.na(g$leftGene))] = 
	g$rightGeneSym[which(is.na(g$newGeneSym) & is.na(g$leftGene))] 
g$newGeneSym[which(is.na(g$newGeneSym) & is.na(g$rightGene))] = 
	g$leftGeneSym[which(is.na(g$newGeneSym) & is.na(g$rightGene))] 
g$newGeneSym[g$newGeneSym==""] = NA
g$newGeneSym[g$newGeneSym=="-"] = NA
anno$newGeneID = g$newGene
anno$newGeneSymbol = g$newGeneSym
anno$isFusion = grepl("-", anno$newGeneID)
anno$newGeneSymbol[anno$code =="InGen"] = anno$Symbol[anno$code =="InGen"]
anno$newGeneID[anno$code =="InGen"] = anno$gencodeGeneID[anno$code =="InGen"]

## extract out jMap
jMap = anno
colnames(mcols(jMap))[which(colnames(mcols(jMap))=="code")] = "Class"
rm(anno)

############ jCounts
jCounts = as.matrix(as.data.frame(juncCounts$countDF))
jCounts = jCounts[names(jMap),metrics$SAMPLE_ID] # ensure lines up
colnames(jCounts) = metrics$SAMPLE_ID  # change from '.' to hyphens if needed

############ jRpkm
bgJ = matrix(rep(colSums(jCounts)), nc = nrow(metrics), 
	nr = nrow(jCounts),	byrow=TRUE)
jRpkm = jCounts/(bgJ/10e6)

rownames(jCounts) = rownames(jRpkm) = names(jMap)
colnames(jRpkm)  = metrics$SAMPLE_ID 
jMap$meanExprs = rowMeans(jRpkm)


# ## sequence of acceptor/donor sites
# left = right = jMap
# end(left) = start(left) +1
# start(right) = end(right) -1
# jMap$leftSeq  = getSeq(Hsapiens, left)
# jMap$rightSeq = getSeq(Hsapiens, right)



### save counts

tosaveCounts = c("metrics", "geneCounts", "geneMap", "exonCounts", "exonMap", "jCounts", "jMap", 
					"txNumReads", "txMap" )
tosaveRpkm = c("metrics", "geneRpkm", "geneMap", "exonRpkm", "exonMap", "jRpkm", "jMap", 
					"txTpm", "txNumReads", "txMap" )

if (exists("erccTPM")) {
	tosaveCounts = c("erccTPM", tosaveCounts)
	tosaveRpkm = c("erccTPM", tosaveRpkm)
}

save(list=ls()[ls() %in% tosaveCounts], compress=TRUE,
	file= file.path(opt$maindir, paste0('rawCounts_', EXPNAME, '_n', N, '.rda')))
save(list=ls()[ls() %in% tosaveRpkm], compress=TRUE,
	file= file.path(opt$maindir, paste0('rpkmCounts_', EXPNAME, '_n', N, '.rda')))

	
## Create RangedSummarizedExperiment objects
rse_jx <- SummarizedExperiment(assays = list('counts' = jCounts),
    rowRanges = jMap, colData = metrics)
save(rse_jx, file = paste0('rse_jx_', EXPNAME, '_n', N, '.Rdata'))
