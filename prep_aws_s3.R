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

## move into subfolders by hand

thecall = "aws s3 cp S3 s3://spatial-dlpfc/ --recursive"
system(thecall)

## make table of links
urls  = "https://spatial-dlpfc.s3.us-east-2.amazonaws.com/"
f = list.files("S3",recur=TRUE)
f = paste0(urls, f)
names(f) = ss(ss(f,"/",5), "_|\\.")

fList = split(f,names(f))
fMat = do.call("rbind", fList)
fDat = as.data.frame(fMat,stringAsFactors=FALSE)
colnames(fDat) = c("h5_filtered", "h5_raw", "image_hi", "image_lo", "loupe")
fDat$SampleID = rownames(fDat)

fDat = fDat[,c(6,1:5)]
write.table(fDat, file = "AWS_File_locations.tsv",
	sep = "\t", row.names=FALSE)