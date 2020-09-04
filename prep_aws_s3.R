###
library(jaffelab)

dir.create("S3")

## mv loupe files
path = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/"
## copy barcodes etc
h5 = list.files(path, pattern = "feature_bc_matrix.h5", recur=TRUE,full=TRUE)
file.copy(h5, "S3/")

hi = list.files(path, pattern = "image.png", recur=TRUE,full=TRUE)
newhi = paste0("S3/", gsub("/","_",ss(hi,"//", 2)))
file.copy(hi,newhi)

## mv hi res images
pheno = read.csv(paste0(path, "20190925_JHU_Lieber_HumanBrain_LP/dlpfc_visium_imageInfo.csv"),
	as.is=TRUE)
pheno$tifPath = paste0(path, "20190925_JHU_Lieber_HumanBrain_LP/", pheno$image_name)
pheno$newTifPath = paste0(path, "../S3/images/", pheno$sample_name, "_full_image.tif")
file.copy(pheno$tifPath,pheno$newTifPath)


## move into subfolders by hand

## put online
thecall = "aws s3 cp S3 s3://spatial-dlpfc/ --recursive"
system(thecall)

#####
## make table of links
urls  = "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/"
f = list.files("S3",recur=TRUE)
f = paste0(urls, f)
names(f) = ss(ss(f,"/",5), "_|\\.")

fList = split(f,names(f))
fMat = do.call("rbind", fList)
fDat = as.data.frame(fMat,stringAsFactors=FALSE)
colnames(fDat) = c("h5_filtered", "h5_raw", "image_full", "image_hi", "image_lo", "loupe")
fDat$SampleID = rownames(fDat)

fDat = fDat[,c(7,1:6)]
write.table(fDat, file = "AWS_File_locations.tsv",
	sep = "\t", row.names=FALSE)