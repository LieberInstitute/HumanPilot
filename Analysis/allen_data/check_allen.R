
library(jaffelab)

## PDF --> Excel using Adobe Acrobat
## Excel --> TSV using Excel 
## scan subsequent text file
x = scan("HBA_ISH_GeneList.txt",what = "character", sep="\n")

## get rows for each table
ind =c(1, grep("Table", x), nrow(x))
names(ind) = gsub("\"", "", gsub(" ", "", ss(x[ind], "\\.", 1)))
names(ind)[c(1,length(ind))] = c("Table1", "Table5")

## get row indicies
indList = split(ind, names(ind))

first = sapply(indList, "[[", 1)
second = c(first[-1]-1, length(x))
indList2 = mapply(function(x,y) x:y, first,second)

## subset
tableList = lapply(indList2, function(ii) x[ii])

## manually clean up
tableList$Table1 = tableList$Table1[-(1:3)]
tableList = lapply(tableList, function(x) x[-(1:2)])

## read in
datList = lapply(tableList[-3], function(x) {
	read.delim(text = x,as.is=TRUE,header=TRUE)
})
datList = lapply(datList, function(x) {
	x[,colMeans(is.na(x)) < 1]
})

## fix table 3
tab3 = tableList$Table3
tab3 = tab3[!grepl("Gene", tab3)]
tab3 = tab3[!grepl("Character", tab3)]
tab3 = tab3[!grepl("Table 3", tab3)]
tab3 = tab3[!grepl("System", tab3)]
tab3 = tab3[!grepl("Symbol", tab3)]

dat3 = read.delim(text = tab3,as.is=TRUE,header=FALSE)
dat3 = dat3[,colMeans(is.na(dat3)) < 1]
colnames(dat3) = c("Gene.Symbol", "EntrezID", "Gene.Description", 
	"Characterized", "System", "Family" , "Marker.Type")
datList$Table3 = dat3

datList = datList[paste0("Table", 1:5)]
dir.create("tables")
for(i in seq(along=datList)) {
	write.table(datList[[i]], paste0("tables/allen_HBA_ISH_GeneList_Table", i, ".txt"),
		row.names=FALSE)
}
zip("allen_HBA_ISH_GeneList_Tables.zip", files=list.files("tables",full=TRUE))

