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

## Re-read this data
# dat <-
#     read.table(here::here("Analysis", "Layer_Guesses", "suppTable_ldsc.txt"),
#         header = TRUE)

options(width = 160)
dat[dat$Coefficient_p < 0.01,]
#    Category              Trait  Prop._SNPs   Prop._h2 Enrichment Enrichment_p Enrichment_holm  Coefficient Coefficient_z.score Coefficient_p Coefficient_holm
# 66   Layer2   Bipolar_disorder 0.006832419 0.05541614   8.110764 1.478401e-07    1.034880e-06 3.107158e-07            4.906768  4.629477e-07     3.240634e-06
# 70   Layer2 College_attainment 0.006832419 0.06533992   9.563220 3.283226e-05    2.298258e-04 9.556632e-08            3.900983  4.790141e-05     3.353099e-04
# 82   Layer2                 IQ 0.006832419 0.03659739   5.356432 8.706680e-03    6.094676e-02 1.148389e-07            2.478732  6.592522e-03     4.614765e-02
# 88   Layer2      Schizophrenia 0.006832419 0.03178507   4.652097 6.418024e-05    4.492617e-04 2.172155e-07            3.605726  1.556408e-04     1.089485e-03

dat[which(dat$Trait == "Schizophrenia"),]
#     Category         Trait  Prop._SNPs     Prop._h2 Enrichment Enrichment_p Enrichment_holm   Coefficient Coefficient_z.score Coefficient_p Coefficient_holm
# 28        WM Schizophrenia 0.023137726  0.035500713  1.5343216 2.809097e-01    0.5618193293 -3.987475e-09          -0.1065047  0.5424090553      1.000000000
# 58    Layer1 Schizophrenia 0.009772554 -0.002301291 -0.2354851 4.182254e-02    0.2091127237 -1.231121e-07          -2.7544016  0.9970600218      1.000000000
# 88    Layer2 Schizophrenia 0.006832419  0.031785074  4.6520969 6.418024e-05    0.0004492617  2.172155e-07           3.6057259  0.0001556408      0.001089485
# 118   Layer3 Schizophrenia 0.001366383  0.002723213  1.9930093 5.126573e-01    0.5618193293  3.616103e-08           0.3473531  0.3641630397      1.000000000
# 148   Layer4 Schizophrenia 0.002765668  0.009106806  3.2928058 1.129455e-01    0.3388364215  1.200315e-07           1.2214042  0.1109665067      0.443866027
# 178   Layer5 Schizophrenia 0.003655723  0.014339281  3.9224200 1.349028e-02    0.0809416942  1.716627e-07           2.1686229  0.0150556616      0.090333970
# 208   Layer6 Schizophrenia 0.002175953  0.009330657  4.2880780 8.047206e-02    0.3218882452  1.932178e-07           1.4978444  0.0670868440      0.335434220

dat[which(dat$Trait == "Schizophrenia" & dat$Category == "Layer2"),]
#    Category         Trait  Prop._SNPs   Prop._h2 Enrichment Enrichment_p Enrichment_holm  Coefficient Coefficient_z.score Coefficient_p Coefficient_holm
# 88   Layer2 Schizophrenia 0.006832419 0.03178507   4.652097 6.418024e-05    0.0004492617 2.172155e-07            3.605726  0.0001556408      0.001089485

library(tidyr)
wide_enr_enrichment = spread(
		dat[,c("Category", "Trait", "Enrichment")],
	Trait, Enrichment)
wide_enr_p = spread(
		dat[,c("Category", "Trait", "Enrichment_p")],
	Trait, Enrichment_p)
rownames(wide_enr_p) = wide_enr_p$Category
wide_enr_p$Category = NULL
options(width = 100)
round(-log10(wide_enr_p),1)
#        ADHD Alzheimers_disease Anorexia_nervosa Anxiety_disorder Autism_spectrum_disorder
# Layer1  0.0                1.2              0.1              0.0                      0.2
# Layer2  1.1                0.8              0.0              0.7                      0.0
# Layer3  0.2                0.3              0.1              0.1                      0.2
# Layer4  0.1                0.3              0.7              1.2                      0.0
# Layer5  0.2                1.0              0.4              0.9                      0.7
# Layer6  0.3                0.2              0.0              2.6                      0.5
# WM      0.5                2.3              0.4              0.3                      1.6
#        Bipolar_disorder BMI Childhood_cognitive_performance Cigarettes_per_day College_attainment
# Layer1              0.3 0.1                             0.3                0.8                0.0
# Layer2              6.8 2.0                             1.0                0.1                4.5
# Layer3              1.6 0.6                             0.6                0.1                0.4
# Layer4              1.1 0.9                             0.1                0.7                0.3
# Layer5              2.0 1.0                             0.2                0.5                0.1
# Layer6              0.8 0.2                             1.0                0.1                0.0
# WM                  0.3 0.8                             0.5                0.2                0.2
#        Conscientiousness Coronary_artery_disease Crohns_disease Depressive_symptoms Epilepsy
# Layer1               0.1                     0.5            0.5                 0.5      1.5
# Layer2               0.9                     1.9            0.6                 0.7      0.8
# Layer3               0.1                     0.6            0.0                 0.6      1.1
# Layer4               0.2                     0.4            1.2                 0.2      1.1
# Layer5               0.4                     0.8            0.6                 1.3      0.1
# Layer6               1.4                     0.3            0.2                 0.3      0.4
# WM                   0.6                     1.2            1.6                 0.0      0.2
#        Ever_smoked Extraversion Focal_epilepsy Generalized_epilepsy Height Intracarebral_hemorrhage
# Layer1         3.0          0.1            1.0                  0.1    3.1                      0.5
# Layer2         0.7          0.0            0.1                  0.2    0.3                      0.9
# Layer3         0.1          0.6            0.3                  0.4    0.0                      0.6
# Layer4         0.2          0.2            0.3                  0.6    0.1                      0.3
# Layer5         0.4          0.5            0.6                  0.6    0.3                      0.0
# Layer6         0.5          0.6            0.4                  0.0    0.9                      1.1
# WM             0.3          0.0            1.1                  0.5    4.4                      1.6
#         IQ Ischemic_stroke Major_depressive_disorder Neuroticism Openness PTSD Schizophrenia
# Layer1 0.1             0.3                       0.3         0.4      0.3  0.0           1.4
# Layer2 2.1             0.1                       0.8         0.2      0.1  0.0           4.2
# Layer3 0.9             0.3                       0.7         0.8      0.2  0.8           0.3
# Layer4 0.2             0.1                       0.4         0.9      0.4  0.2           0.9
# Layer5 0.7             0.6                       1.6         1.1      0.1  0.3           1.9
# Layer6 0.4             1.0                       0.9         0.8      0.3  0.4           1.1
# WM     0.6             0.2                       0.9         0.1      0.2  0.0           0.6
#        Subjective_well-being Years_of_education
# Layer1                   0.4                0.1
# Layer2                   0.4                1.0
# Layer3                   0.2                0.4
# Layer4                   0.5                0.1
# Layer5                   0.7                0.1
# Layer6                   0.2                0.1
# WM                       0.1                0.2
round(-log10(wide_enr_p),1) > -log10(0.05)
#         ADHD Alzheimers_disease Anorexia_nervosa Anxiety_disorder Autism_spectrum_disorder Bipolar_disorder   BMI
# Layer1 FALSE              FALSE            FALSE            FALSE                    FALSE            FALSE FALSE
# Layer2 FALSE              FALSE            FALSE            FALSE                    FALSE             TRUE  TRUE
# Layer3 FALSE              FALSE            FALSE            FALSE                    FALSE             TRUE FALSE
# Layer4 FALSE              FALSE            FALSE            FALSE                    FALSE            FALSE FALSE
# Layer5 FALSE              FALSE            FALSE            FALSE                    FALSE             TRUE FALSE
# Layer6 FALSE              FALSE            FALSE             TRUE                    FALSE            FALSE FALSE
# WM     FALSE               TRUE            FALSE            FALSE                     TRUE            FALSE FALSE
#        Childhood_cognitive_performance Cigarettes_per_day College_attainment Conscientiousness Coronary_artery_disease
# Layer1                           FALSE              FALSE              FALSE             FALSE                   FALSE
# Layer2                           FALSE              FALSE               TRUE             FALSE                    TRUE
# Layer3                           FALSE              FALSE              FALSE             FALSE                   FALSE
# Layer4                           FALSE              FALSE              FALSE             FALSE                   FALSE
# Layer5                           FALSE              FALSE              FALSE             FALSE                   FALSE
# Layer6                           FALSE              FALSE              FALSE              TRUE                   FALSE
# WM                               FALSE              FALSE              FALSE             FALSE                   FALSE
#        Crohns_disease Depressive_symptoms Epilepsy Ever_smoked Extraversion Focal_epilepsy Generalized_epilepsy Height
# Layer1          FALSE               FALSE     TRUE        TRUE        FALSE          FALSE                FALSE   TRUE
# Layer2          FALSE               FALSE    FALSE       FALSE        FALSE          FALSE                FALSE  FALSE
# Layer3          FALSE               FALSE    FALSE       FALSE        FALSE          FALSE                FALSE  FALSE
# Layer4          FALSE               FALSE    FALSE       FALSE        FALSE          FALSE                FALSE  FALSE
# Layer5          FALSE               FALSE    FALSE       FALSE        FALSE          FALSE                FALSE  FALSE
# Layer6          FALSE               FALSE    FALSE       FALSE        FALSE          FALSE                FALSE  FALSE
# WM               TRUE               FALSE    FALSE       FALSE        FALSE          FALSE                FALSE   TRUE
#        Intracarebral_hemorrhage    IQ Ischemic_stroke Major_depressive_disorder Neuroticism Openness  PTSD Schizophrenia
# Layer1                    FALSE FALSE           FALSE                     FALSE       FALSE    FALSE FALSE          TRUE
# Layer2                    FALSE  TRUE           FALSE                     FALSE       FALSE    FALSE FALSE          TRUE
# Layer3                    FALSE FALSE           FALSE                     FALSE       FALSE    FALSE FALSE         FALSE
# Layer4                    FALSE FALSE           FALSE                     FALSE       FALSE    FALSE FALSE         FALSE
# Layer5                    FALSE FALSE           FALSE                      TRUE       FALSE    FALSE FALSE          TRUE
# Layer6                    FALSE FALSE           FALSE                     FALSE       FALSE    FALSE FALSE         FALSE
# WM                         TRUE FALSE           FALSE                     FALSE       FALSE    FALSE FALSE         FALSE
#        Subjective_well-being Years_of_education
# Layer1                 FALSE              FALSE
# Layer2                 FALSE              FALSE
# Layer3                 FALSE              FALSE
# Layer4                 FALSE              FALSE
# Layer5                 FALSE              FALSE
# Layer6                 FALSE              FALSE
# WM                     FALSE              FALSE

wide_enr_cp = spread(
		dat[,c("Category", "Trait", "Coefficient_p")],
	Trait, Coefficient_p)
rownames(wide_enr_cp) = wide_enr_cp$Category
wide_enr_cp$Category = NULL
round(-log10(wide_enr_cp),1)
#        ADHD Alzheimers_disease Anorexia_nervosa Anxiety_disorder Autism_spectrum_disorder Bipolar_disorder BMI
# Layer1  0.3                1.0              0.4              0.2                      0.5              0.2 0.0
# Layer2  1.4                0.5              0.3              0.9                      0.2              6.3 1.9
# Layer3  0.5                0.1              0.4              0.2                      0.1              1.7 0.0
# Layer4  0.4                0.4              1.0              1.5                      0.3              1.0 0.8
# Layer5  0.5                1.1              0.7              1.2                      0.0              1.8 0.9
# Layer6  0.1                0.1              0.3              0.0                      0.8              0.8 0.3
# WM      1.1                1.1              0.8              0.1                      1.7              0.0 0.5
#        Childhood_cognitive_performance Cigarettes_per_day College_attainment Conscientiousness Coronary_artery_disease
# Layer1                             0.1                0.0                0.2               0.3                     0.5
# Layer2                             1.2                0.3                4.3               1.1                     1.8
# Layer3                             0.9                0.2                0.6               0.2                     0.8
# Layer4                             0.4                0.0                0.5               0.4                     0.6
# Layer5                             0.1                0.1                0.3               0.6                     0.9
# Layer6                             1.3                0.4                0.3               0.0                     0.5
# WM                                 0.7                0.2                0.1               0.6                     0.7
#        Crohns_disease Depressive_symptoms Epilepsy Ever_smoked Extraversion Focal_epilepsy Generalized_epilepsy Height
# Layer1            0.1                 0.7      0.0         0.0          0.4            0.0                  0.1    1.7
# Layer2            0.4                 0.8      0.9         0.0          0.3            0.2                  0.4    0.1
# Layer3            0.2                 0.8      1.4         0.2          0.1            0.6                  0.7    0.1
# Layer4            0.9                 0.4      1.3         0.5          0.2            0.6                  0.8    0.1
# Layer5            0.5                 1.4      0.2         0.7          0.1            0.1                  0.0    0.0
# Layer6            0.3                 0.6      0.7         0.9          0.0            0.7                  0.3    0.5
# WM                0.4                 0.1      0.1         0.8          0.3            0.0                  0.7    1.3
#        Intracarebral_hemorrhage  IQ Ischemic_stroke Major_depressive_disorder Neuroticism Openness PTSD Schizophrenia
# Layer1                      0.7 0.2             0.5                       0.0         0.1      0.1  0.5           0.0
# Layer2                      1.1 2.2             0.3                       0.7         0.5      0.3  0.4           3.8
# Layer3                      0.9 1.1             0.1                       0.9         1.0      0.2  1.2           0.4
# Layer4                      0.6 0.1             0.4                       0.5         1.2      0.8  0.2           1.0
# Layer5                      0.3 0.9             0.0                       1.6         1.4      0.4  0.7           1.8
# Layer6                      0.0 0.6             0.0                       1.0         1.1      0.6  0.1           1.2
# WM                          1.5 0.7             0.3                       0.5         0.3      0.2  0.4           0.3
#        Subjective_well-being Years_of_education
# Layer1                   0.1                0.3
# Layer2                   0.5                1.0
# Layer3                   0.4                0.5
# Layer4                   0.1                0.2
# Layer5                   0.9                0.2
# Layer6                   0.2                0.3
# WM                       0.1                0.1

round(-log10(wide_enr_cp),1) > -log10(0.05)
#         ADHD Alzheimers_disease Anorexia_nervosa Anxiety_disorder Autism_spectrum_disorder Bipolar_disorder   BMI
# Layer1 FALSE              FALSE            FALSE            FALSE                    FALSE            FALSE FALSE
# Layer2  TRUE              FALSE            FALSE            FALSE                    FALSE             TRUE  TRUE
# Layer3 FALSE              FALSE            FALSE            FALSE                    FALSE             TRUE FALSE
# Layer4 FALSE              FALSE            FALSE             TRUE                    FALSE            FALSE FALSE
# Layer5 FALSE              FALSE            FALSE            FALSE                    FALSE             TRUE FALSE
# Layer6 FALSE              FALSE            FALSE            FALSE                    FALSE            FALSE FALSE
# WM     FALSE              FALSE            FALSE            FALSE                     TRUE            FALSE FALSE
#        Childhood_cognitive_performance Cigarettes_per_day College_attainment Conscientiousness Coronary_artery_disease
# Layer1                           FALSE              FALSE              FALSE             FALSE                   FALSE
# Layer2                           FALSE              FALSE               TRUE             FALSE                    TRUE
# Layer3                           FALSE              FALSE              FALSE             FALSE                   FALSE
# Layer4                           FALSE              FALSE              FALSE             FALSE                   FALSE
# Layer5                           FALSE              FALSE              FALSE             FALSE                   FALSE
# Layer6                           FALSE              FALSE              FALSE             FALSE                   FALSE
# WM                               FALSE              FALSE              FALSE             FALSE                   FALSE
#        Crohns_disease Depressive_symptoms Epilepsy Ever_smoked Extraversion Focal_epilepsy Generalized_epilepsy Height
# Layer1          FALSE               FALSE    FALSE       FALSE        FALSE          FALSE                FALSE   TRUE
# Layer2          FALSE               FALSE    FALSE       FALSE        FALSE          FALSE                FALSE  FALSE
# Layer3          FALSE               FALSE     TRUE       FALSE        FALSE          FALSE                FALSE  FALSE
# Layer4          FALSE               FALSE    FALSE       FALSE        FALSE          FALSE                FALSE  FALSE
# Layer5          FALSE                TRUE    FALSE       FALSE        FALSE          FALSE                FALSE  FALSE
# Layer6          FALSE               FALSE    FALSE       FALSE        FALSE          FALSE                FALSE  FALSE
# WM              FALSE               FALSE    FALSE       FALSE        FALSE          FALSE                FALSE  FALSE
#        Intracarebral_hemorrhage    IQ Ischemic_stroke Major_depressive_disorder Neuroticism Openness  PTSD Schizophrenia
# Layer1                    FALSE FALSE           FALSE                     FALSE       FALSE    FALSE FALSE         FALSE
# Layer2                    FALSE  TRUE           FALSE                     FALSE       FALSE    FALSE FALSE          TRUE
# Layer3                    FALSE FALSE           FALSE                     FALSE       FALSE    FALSE FALSE         FALSE
# Layer4                    FALSE FALSE           FALSE                     FALSE       FALSE    FALSE FALSE         FALSE
# Layer5                    FALSE FALSE           FALSE                      TRUE        TRUE    FALSE FALSE          TRUE
# Layer6                    FALSE FALSE           FALSE                     FALSE       FALSE    FALSE FALSE         FALSE
# WM                         TRUE FALSE           FALSE                     FALSE       FALSE    FALSE FALSE         FALSE
#        Subjective_well-being Years_of_education
# Layer1                 FALSE              FALSE
# Layer2                 FALSE              FALSE
# Layer3                 FALSE              FALSE
# Layer4                 FALSE              FALSE
# Layer5                 FALSE              FALSE
# Layer6                 FALSE              FALSE
# WM                     FALSE              FALSE
