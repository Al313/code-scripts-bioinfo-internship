# In this script I want to count the number of snvs for each sample

sample_dir <- "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF/annotated-variants/cancer-type/"


metadata <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/labelled_metadata2.tsv", sep = "\t", header = T, stringsAsFactors = F)



metadata$snv_count <- NA
for (i in 1:nrow(metadata)){
  message(i)
  vcf <- read.csv(file = paste0(sample_dir, metadata$primaryTumorLocation[i], "/", metadata$sampleId[i], "/annotated-", metadata$sampleId[i], ".txt"), sep = "\t", header = T, stringsAsFactors = F)
  metadata$snv_count[i] <- nrow(vcf)

}


saveRDS(metadata, "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/metadata_with_snv_counts.rds")









# 
# metadata <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/metadata_with_snv_counts.rds")
# quantile(metadata$snv_count, probs = .9625)
# mean(metadata$snv_count)/3000
# median(metadata$snv_count)/3000
# 
# boxplot(log10(metadata$snv_count), pch = 19, main = "Distribution of SNV counts in HMF dataset", ylab = "Log10(SNV counts)")
# boxplot(log10(metadata$snv_count[metadata$pathway == "Grey-listed"]), pch = 19, main = "Distribution of SNV counts in HMF dataset", ylab = "Log10(SNV counts)")
# 


