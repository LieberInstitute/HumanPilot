# Plot LDSC results for 'adjusting for baseline' analyses
# code adapted from Peter Hickey (https://github.com/hansenlab/BrainEpigenomeNN/blob/master/SLDSR/scripts/plot_ldsc.baseline_adjustments.R)
# qrsh -l bluejay,mem_free=100G,h_vmem=100G,h_fsize=500G

### ----------------------------------------------------------------------------
### Packages
###

library(GenomicRanges)
library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(gplots)
library(cowplot)
library(RColorBrewer)
library(jaffelab)
library(parallel)

###-----------------------------------------------------------------------------
### Category sizes
###

# ------------------------------------------------------------------------------
# Brain categories
#

brain_categories <- readRDS("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/LDSC/categories.rds")

brain_categories_df <- data.frame(Category = names(brain_categories),stringsAsFactors=FALSE)
brain_categories_df$Layer = ss(brain_categories_df$Category, "_")
brain_categories_df$Type = ifelse(grepl("Neuropil", brain_categories_df$Category), "Neuropil", "Laminar")


## ----------------------------------------------------------------------------
### Load data and construct objects
###

fls <- unlist(lapply(names(brain_categories), function(bc) {
  fls <- list.files("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/LDSC/output",
                    pattern = glob2rx(
                      paste0(bc, ".*Phase1.results")),
                    full.names = TRUE)
  grep("adjusting", fls, invert = TRUE, value = TRUE)
}))

## drop neuropil
fls = fls[!grepl("Neuropil", fls)]

# Read in files, tidy up, and rbind
x <- bind_rows(mclapply(fls, function(fl) {
  cat(".")
  suppressMessages(read_tsv(fl)) %>%
    filter(Category == "L2_0" | Category == "CNS_0") %>%
    mutate(Category = sapply(strsplit(basename(fl), "\\."), "[[", 1),
           Trait = sapply(strsplit(sub(".Phase1.results", "", basename(fl)),
                                   "\\."),
                          "[[", 2),
           lower = Enrichment - 1.96 * Enrichment_std_error,
           upper = Enrichment + 1.96 * Enrichment_std_error,
           file = fl) %>%
    mutate(Category = factor(Category, Category, ordered = TRUE))
}, mc.cores=10))

# NOTE: Anttila report these traits "had in sufficient evidence of additive
#       heritability for robust analysis" and excluded them from further
#       analysis
x <- x %>%
  filter(Trait != "Agreeableness",
         Trait != "Cardioembolic_stroke",
         Trait != "Large-vessel_disease",
         Trait != "Small-vessel_disease")
fls <- fls[-c(grep("Agreeableness", fls), grep("Cardioembolic_stroke", fls), 
              grep("Large-vessel_disease", fls), grep("Small-vessel_disease", fls))]

# also use updated summary statistics for 3 psychiatric diseases since Feinberg paper
x <- x %>%
  filter(Trait != "Autism_spectrum_disorder",
         Trait != "Bipolar_disorder",
         Trait != "Major_depressive_disorder")
fls <- fls[-c(grep("Autism_spectrum_disorder.", fls, fixed=T), grep("Bipolar_disorder.", fls, fixed=T), 
              grep("Major_depressive_disorder.", fls, fixed=T))]


# Join munged LDSC output with categories_df and traits_df

x$Trait = gsub("_latest","", x$Trait)

x <- x %>%
  group_by(Category, file) %>%
  filter(grepl(Category, file)) %>%
  ungroup() %>%
  mutate(Coefficient_p = pnorm(`Coefficient_z-score`, lower.tail = FALSE)) %>%
  group_by(Trait) %>%
  mutate(Coefficient_holm = p.adjust(Coefficient_p, method = "holm"),
         Enrichment_holm = p.adjust(Enrichment_p, method = "holm")) %>%
  ungroup()

dat = as.data.frame(x)
dat = dat[,c(1,11,2:3,5,7,17,8,10,15:16)]
write.table(dat, file = "../suppTable_ldsc.txt",row.names=FALSE,sep="\t")

dat[dat$Coefficient_p < 0.01,]

dat[which(dat$Trait == "Schizophrenia")[1:13],]
dat[which(dat$Trait == "Schizophrenia" & dat$Category == "Layer2"),]

library(tidyr)
wide_enr_enrichment = spread(
		dat[,c("Category", "Trait", "Enrichment")],
	Trait, Enrichment)
wide_enr_p = spread(
		dat[,c("Category", "Trait", "Enrichment_p")],
	Trait, Enrichment_p)
rownames(wide_enr_p) = wide_enr_p$Category
wide_enr_p$Category = NULL
round(-log10(wide_enr_p),1)

wide_enr_cp = spread(
		dat[,c("Category", "Trait", "Coefficient_p")],
	Trait, Coefficient_p)
rownames(wide_enr_cp) = wide_enr_cp$Category
wide_enr_cp$Category = NULL
round(-log10(wide_enr_cp),1)
