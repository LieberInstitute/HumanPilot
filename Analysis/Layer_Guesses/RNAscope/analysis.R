##
library(RColorBrewer)
library(spatialLIBD)

## read in data
man = read.csv("/dcl01/lieber/ajaffe/Maddy/RNAscope/20x_Kristen_Visium_paper/Subject3_AQP4570_RELN520_TRABD2A620_BCL11B690_20xtile_Linear unmixing_Stitch_Image.csv",
	as.is=TRUE)
dat = read.csv("/dcl01/lieber/ajaffe/Maddy/RNAscope/20x_Kristen_Visium_paper/Subject3_AQP4570_RELN520_TRABD2A620_BCL11B690_20xtile_Linear unmixing_Stitch_ROI.csv",
	as.is=TRUE)

## add proportion of pixels
dat$PP_O69_Opal520_20x = dat$MP_O69_Opal520_20x / dat$RVolume
dat$PP_O68_Opal570_20x = dat$MP_O68_Opal570_20x / dat$RVolume
dat$PP_O70_Opal620_20x = dat$MP_O70_Opal620_20x / dat$RVolume
dat$PP_O71_Opal690_20x = dat$MP_O71_Opal690_20x / dat$RVolume

plot(dat$CentroidX, dat$CentroidY)

## manually change per image
dat$CentroidY_new= dat$CentroidY-0.4*dat$CentroidX
plot(dat$CentroidX, dat$CentroidY_new)

dat$Layer_Group = cut(dat$CentroidY_new, 
	c(0,2000, 5000, 6500,8000, 10000,13000),
	labels = paste0("Layer", 6:1))
dat$Layer_Group = factor(as.character(dat$Layer_Group),
	paste0("Layer", 1:6))
	
dat$Layer_Colors = libd_layer_colors[as.character(dat$Layer_Group)]

palette(brewer.pal(7,"Dark2"))
plot(dat$CentroidX, dat$CentroidY_new, 
	pch = 19, col = dat$Layer_Colors )	
## RELN: 520
## AQP4: 570
## TRABD2A: 620
## BCL11B: 690

## RELN
table( dat$Layer_Group, dat$MD_O69_Opal520_20x > 2)
prop.table(table( dat$Layer_Group, 
		dat$MD_O69_Opal520_20x > 4), 1)
prop.table(table( dat$Layer_Group, 
		dat$PP_O69_Opal520_20x > 0.05), 1)

## AQP4
prop.table(table( dat$Layer_Group, 
		dat$MD_O68_Opal570_20x > 4), 1)
prop.table(table( dat$Layer_Group, 
		dat$PP_O68_Opal570_20x > 0.05), 1)

## TRABD2A
prop.table(table( dat$Layer_Group, 
		dat$MD_O70_Opal620_20x > 4), 1)
prop.table(table( dat$Layer_Group, 
		dat$PP_O70_Opal620_20x > 0.05), 1)


palette(c("white",brewer.pal(5, "Blues")))
plot(dat$CentroidX, dat$CentroidY_new,pch=19,
	col = floor(sqrt(dat$MD_O69_Opal520_20x))+1)
plot(dat$CentroidX, dat$CentroidY_new,pch=19,
	col = floor(sqrt(dat$MD_O68_Opal570_20x))+1)
plot(dat$CentroidX, dat$CentroidY_new,pch=19,
	col = floor(sqrt(dat$MD_O70_Opal620_20x))+1)
plot(dat$CentroidX, dat$CentroidY_new,pch=19,
	col = floor(sqrt(dat$MD_O71_Opal690_20x))+1)

plot(MD_O69_Opal520_20x ~ CentroidY, data = dat)
plot(PP_O68_Opal570_20x ~ cut(CentroidY_new,7), data = dat)

plot(MD_O70_Opal620_20x ~ cut(CentroidY,7), data = dat)
plot(PP_O70_Opal620_20x ~ cut(CentroidY,7), data = dat)

plot(MD_O71_Opal690_20x ~ cut(CentroidY,7), data = dat)
plot(PP_O71_Opal690_20x ~ cut(CentroidY,7), data = dat)

