options(width=100)
library(scater)
library(SingleCellExperiment)

## load annotated sce objects
load("SCE_singlet-spots_MNT.rda")
load("SCE_neuropil-spots_MNT.rda")

dim(sce.singlet)
dim(sce.neuropil)

table(sce.singlet$subject_position, sce.singlet$prelimCluster_MNT)
table(sce.neuropil$subject_position, sce.neuropil$prelimCluster_MNT)

## load layer guesses
load("../Layer_Guesses/rda/layer_guess_tab.Rdata")

## add layer info
layer_guess_tab$layer[layer_guess_tab$layer == 'Layer 2/3'] <- 'Layer 3'
layer_guess_tab$Layer = gsub("ayer ", "", layer_guess_tab$layer)

## add
sce.singlet$Layer = layer_guess_tab$Layer[match(sce.singlet$key, layer_guess_tab$key)]
sce.neuropil$Layer = layer_guess_tab$Layer[match(sce.neuropil$key, layer_guess_tab$key)]

## cross tables
tt.singlet = table(sce.singlet$prelimCluster_MNT, sce.singlet$Layer)
tt.singlet
round(prop.table(tt.singlet,1)*100,1)

tt.neuropil = table(sce.neuropil$prelimCluster_MNT, sce.neuropil$Layer)
round(prop.table(tt.neuropil,1)*100,1)

table(sce.neuropil$prelimCluster_MNT, sce.neuropil$Layer)