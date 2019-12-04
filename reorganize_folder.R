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
