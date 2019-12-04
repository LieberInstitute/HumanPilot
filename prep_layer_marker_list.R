###
###
library(readxl)
library(biomaRt)

## read in genes
geneTab = read_excel("cortical layer marker gene list_1.xlsx")
geneTab = as.data.frame(geneTab)
colnames(geneTab)[1] = "Gene"

## get human/mouse conversion
ensembl=useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
MMtoHG = getBM(attributes = c('ensembl_gene_id','external_gene_name', 'hsapiens_homolog_ensembl_gene'),mart = ensembl)

## fix using ensembl
geneTab$MouseEnsID = MMtoHG$ensembl_gene_id[match(tolower(geneTab$Gene), tolower(MMtoHG$external_gene_name))]

## fix up
geneTab$MouseEnsID[geneTab$Gene == "A930038C07Rik"] = "ENSMUSG00000049001"
geneTab$MouseEnsID[geneTab$Gene == "9830123M21Rik"] = "ENSMUSG00000051985"
geneTab$MouseEnsID[geneTab$Gene == "Mirn189"] = "ENSMUSG00000105904"
geneTab$MouseEnsID[geneTab$Gene == "kitlg"] = "ENSMUSG00000019966"
geneTab$MouseEnsID[geneTab$Gene == "Brn1"] = "ENSMUSG00000045515"
geneTab$MouseEnsID[geneTab$Gene == "Brn2"] = "ENSMUSG00000095139"
geneTab$MouseEnsID[geneTab$Gene == "C030003D03Rik"] = "ENSMUSG00000049420"
geneTab$MouseEnsID[geneTab$Gene == "Trb"] = "ENSMUSG00000018697"
geneTab$MouseEnsID[geneTab$Gene == "Igh6"] = "ENSMUSG00000076617"

## rematch
geneTab$MouseSym = MMtoHG$external_gene_name[match(geneTab$MouseEnsID, MMtoHG$ensembl_gene_id)]
geneTab$HumanEnsID = MMtoHG$hsapiens_homolog_ensembl_gene[match(geneTab$MouseEnsID, MMtoHG$ensembl_gene_id)]

## layer labels
geneTab$Label = apply(geneTab[,2:9], 1, function(x) paste(colnames(geneTab)[2:9][which(x!=0)],collapse="/"))
geneTab$MouseSym[geneTab$MouseSym == "6430573F11Rik"] = "643..Rik"
geneTab = geneTab[!is.na(geneTab$HumanEnsID),]
write.csv(geneTab, file = "mouse_layer_marker_info_cleaned.csv",row.names=FALSE)