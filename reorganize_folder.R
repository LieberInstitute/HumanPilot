###
library(jaffelab)
library(filesstrings)

## unpack
system("tar xvf Lieber_Visium_Transfer_1.tar")
system("tar xvf Lieber_Visium_Transfer_2.tar")

## move to common folder
system("mv Lieber_Visium_Transfer_1/ 10X")
system("mv Lieber_Visium_Transfer_2/ 10X")
system("mv Lieber_Visium_Transfer_2/Lieber/ Analysis/")

## make sample-specific folders
path = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X"
f = list.files(path, full =TRUE)
names(f) = list.files(path)

fList = split(f, ss(names(f), "_"))

for(i in seq(along=fList)) {
	newpath = paste0(path, "/", names(fList)[i])
	dir.create(newpath)
	file.move(fList[[i]], newpath) 
}

## move analysis output
path2 = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis"
f2 = list.files(path2,full =TRUE,recur=TRUE)
names(f2) = ss(f2, "/", 8)
f2 = f2[grepl("^1", names(f2))]

fList2 = split(f2, names(f2))

for(i in seq(along=fList2)) {
	newpath = paste0(path, "/", names(fList)[i])
	file.move(fList2[[i]], newpath) 
}

### and loupe files
path3 = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/Lieber_websummary_and_loupe"
f3 = list.files(path3,full =TRUE,recur=TRUE)
names(f3) = ss(ss(f3, "/", 9), "\\.")
names(f3) = ss(names(f3), "_")

fList3 = split(f3, names(f3))

for(i in seq(along=fList3)) {
	newpath = paste0(path, "/", names(fList3)[i])
	file.move(fList3[[i]], newpath,overwrite=TRUE) 
}