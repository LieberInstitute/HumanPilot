###
library(jaffelab)
metricFiles = list.files("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X",
	pattern = "metrics_summary_csv.csv",full=TRUE,recur=TRUE)
names(metricFiles) = ss(metricFiles, "/", 8)

metrics = sapply(metricFiles, read.csv,as.is=TRUE)
write.table(metrics,sep = "\t",
	file = "visium_dlpfc_pilot_sample_metrics.tsv")
