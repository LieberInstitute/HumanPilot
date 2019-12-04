####

library(EBImage)

tifs = paste0("/dcl01/lieber/ajaffe/Maddy/test", 1:6, ".tiff")
x = readImage(tifs[1], info = TRUE)
display(x, method = "raster", all = TRUE)

## make grey
xGrey = x
colorMode(xGrey) = Grayscale
display(xGrey, method = "raster", all = TRUE)
