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
manifest <- read.table('/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/10X/sample_level_layer_map.tsv', 
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
geneFn <- paste0("/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/10X/",
	manifest$SlideID, "/Layers/", manifest$SampleID, ".genes.counts")
names(geneFn) = manifest$SampleID
stopifnot(all(file.exists(geneFn)))

### read in annotation ##
geneMap = read.delim(geneFn[1], skip=1, as.is=TRUE)[,1:6]

## organize gene map
geneMap$Chr = ss(geneMap$Chr, ";")
geneMap$Chr = paste0("chr", geneMap$Chr)
geneMap$Start = as.numeric(ss(geneMap$Start, ";"))
tmp = strsplit(geneMap$End, ";")
geneMap$End = as.numeric(sapply(tmp, function(x) x[length(x)]))
geneMap$Strand = ss(geneMap$Strand, ";")
rownames(geneMap) = geneMap$Geneid
geneMap$gencodeID = geneMap$Geneid
geneMap$ensemblID = ss(geneMap$Geneid, "\\.")
geneMap$Geneid = NULL
geneMap$gene_type = gencodeGENES[geneMap$gencodeID,"gene_biotype"]
geneMap$Symbol = gencodeGENES[geneMap$gencodeID,"gene_name"]


## counts
geneCountList = mclapply(geneFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=4)
geneCounts = do.call("cbind", geneCountList)
rownames(geneCounts) = rownames(geneMap)
geneCounts = geneCounts[,manifest$SampleID] # put in order

# check
colSums(geneCounts)/1e6

# number of reads assigned
geneStatList = mclapply(paste0(geneFn, ".summary"), 
	read.delim,row.names=1, mc.cores=4)
geneStats = do.call("cbind", geneStatList)
colnames(geneStats) = manifest$SampleID
## Add all the other stats from featureCounts at the gene level
geneStats_t <- t(geneStats)
colnames(geneStats_t) <- paste0('gene_', colnames(geneStats_t))
manifest <- cbind(manifest, geneStats_t[,1:4])
manifest$totalAssignedGene = as.numeric(geneStats[1,] / colSums(geneStats))

## Create gene,exon RangedSummarizedExperiment objects
gr_genes <- GRanges(seqnames = geneMap$Chr,
    IRanges(geneMap$Start, geneMap$End), strand = geneMap$Strand)
names(gr_genes) <- rownames(geneMap)
mcols(gr_genes) <- DataFrame(geneMap[, - which(colnames(geneMap) %in%
    c('Chr', 'Start', 'End', 'Strand'))])

rse_gene <- SummarizedExperiment(assays = list('counts' = geneCounts),
    rowRanges = gr_genes, colData = manifest)
save(rse_gene, file = 'rse_gene_layerLevel_n76.Rdata')


###############
### exon counts
exonFn <- paste0("/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/10X/",
	manifest$SlideID, "/Layers/", manifest$SampleID, ".exons.counts")
names(exonFn) = manifest$SampleID
stopifnot(all(file.exists(exonFn)))

### read in annotation ##
exonMap = read.delim(exonFn[1], skip=1, as.is=TRUE)[,1:6]
exonMap$gencodeID = exonMap$Geneid
exonMap$ensemblID = ss(exonMap$Geneid, "\\.")
rownames(exonMap) = paste0("e", rownames(exonMap))
exonMap$Geneid = NULL
exonMap$gene_type = gencodeGENES[exonMap$gencodeID,"gene_biotype"]
exonMap$Symbol = gencodeGENES[exonMap$gencodeID,"gene_name"]
## add gencode exon id
exonMap = join(exonMap, gencodeEXONS, type="left", match="first")
exonMap$Chr = paste0("chr", exonMap$Chr)

## counts
exonCountList = mclapply(exonFn, function(x) {
	cat(".")
	read.delim(pipe(paste("cut -f7", x)), as.is=TRUE,skip=1)[,1]
}, mc.cores=4)
exonCounts = do.call("cbind", exonCountList)
rownames(exonCounts) = rownames(exonMap)
exonCounts = exonCounts[,manifest$SampleID] # put in order

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

gr_exons <- GRanges(seqnames = exonMap$Chr,
    IRanges(exonMap$Start, exonMap$End), strand = exonMap$Strand)
names(gr_exons) <- rownames(exonMap)
mcols(gr_exons) <- DataFrame(exonMap[, - which(colnames(exonMap) %in%
    c('Chr', 'Start', 'End', 'Strand'))])
    
rse_exon <- SummarizedExperiment(assays = list('counts' = exonCounts),
    rowRanges = gr_exons, colData = manifest)
save(rse_exon, file = 'rse_exon_layerLevel_n76.Rdata')

###################
##### junctions

junctionCount = function (junctionFiles, sampleNames = names(junctionFiles),
    output = c("Count", "Rail"), minOverhang = 0, strandSpecific = FALSE,
    illuminaStranded = FALSE, minCount = 1, maxCores = 1) {
    
	theData <- mclapply(junctionFiles, function(x) {
       y <- read.delim(x, skip = 1, header = FALSE,
		  col.names = c("chr", "start", "end", "strand",
			"count"), colClasses = c("character", "integer",
			"integer", "character", "integer"))
        y <- y[y$count >= minCount, ]
        weird <- which(y$strand == "?")
        if (length(weird) > 0) y$strand[weird] = "*"
        gr <- GRanges(y$chr, IRanges(y$start, y$end), strand = y$strand,
					count = y$count)
        return(gr)
        }, mc.cores = maxCores)
  
    
    message(paste(Sys.time(), "creating master table of junctions"))
    grList <- GRangesList(theData)
    if (illuminaStranded & strandSpecific) {
        grList <- GRangesList(mclapply(grList, function(x) {
            strand(x) = ifelse(strand(x) == "+", "-", "+")
            return(x)
        }, mc.cores = maxCores))
    }
    fullGR <- unlist(grList)
    if (!strandSpecific)
        strand(fullGR) <- "*"
    fullGR <- fullGR[!duplicated(fullGR)]
    fullGR <- sort(fullGR)
    fullGR$count <- NULL
    message(paste(Sys.time(), "there are", length(fullGR), "total junctions"))
    message(paste(Sys.time(), "populating count matrix"))
    jNames <- paste0(as.character(seqnames(fullGR)), ":", start(fullGR),
        "-", end(fullGR), "(", as.character(strand(fullGR)),
        ")")
    options(warn = -1)
    mList <- mclapply(grList, match, fullGR, ignore.strand = !strandSpecific,
        mc.cores = maxCores)
    options(warn = 0)
    countList <- mList
    M <- length(jNames)
    message(paste(Sys.time(), "filling in the count matrix"))
    for (i in seq(along = grList)) {
        if (i%%25 == 0)
            cat(".")
        cc <- rep(0, M)
        cc[mList[[i]]] <- theData[[i]]$count
        countList[[i]] <- Rle(cc)
    }
    countDF <- DataFrame(countList, row.names = jNames, check.names = FALSE)
    names(fullGR) <- jNames
    out <- list(countDF = countDF, anno = fullGR)
    return(out)
}


## via primary alignments only
junctionFiles <- paste0("/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/10X/",
	manifest$SlideID, "/Layers/", manifest$SampleID, ".junctions.count")
names(junctionFiles) = manifest$SampleID
	k= file.exists(junctionFiles)

## Read in counts
juncCounts = junctionCount(junctionFiles, manifest$SampleID,
	output = "Count", maxCores=4,strandSpecific=FALSE)

############ anno/jMap
anno = juncCounts$anno
# seqlevels(anno, force=TRUE) = paste0("chr", c(1:22,"X","Y","M"))
seqlevels(anno, pruning.mode="coarse") = c(1:22,"X","Y","MT") # for updated R 3.
seqlevels(anno) = paste0("chr", seqlevels(anno))

## load annotation
load("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/junction_txdb/junction_annotation_hg38_cellranger.rda")
seqlevels(theJunctions)[25] = "chrMT"

## add additional annotation
anno$inGencode = countOverlaps(anno, theJunctions, type="equal") > 0
anno$inGencodeStart = countOverlaps(anno, theJunctions, type="start") > 0
anno$inGencodeEnd = countOverlaps(anno, theJunctions, type="end") > 0

oo = findOverlaps(anno, theJunctions, type="equal")
anno$gencodeGeneID = NA
anno$gencodeGeneID[queryHits(oo)] = as.character(theJunctions$cellrangerID[subjectHits(oo)])
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

############ jCounts
jCounts = as.matrix(as.data.frame(juncCounts$countDF))
jCounts = jCounts[names(jMap),match(manifest$SampleID, colnames(jCounts))] # ensure lines up
colnames(jCounts) = manifest$SampleID  # change from '.' to hyphens if needed

rownames(jCounts) = names(jMap) =  paste0("chr", names(jMap))
	
## Create RangedSummarizedExperiment objects
rse_jx <- SummarizedExperiment(assays = list('counts' = jCounts),
    rowRanges = jMap, colData = manifest)
	
## light filter
keep_jxn = rowSums(assays(rse_jx)$counts) > 1
save(rse_jx, file = 'rse_jx_layerLevel_n76.Rdata')
