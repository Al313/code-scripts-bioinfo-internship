sample_ids <- scan("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/sample-ids/annotated-lung-ids.txt", what = character())

sample_id <- sample_ids[1]

vcf <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF/annotated-variants/lung/",
                              sample_id, "/annotated-", sample_id, ".txt"),
                header = T, sep = "\t", stringsAsFactors = F)
head(vcf)


table(vcf$genomic_func)

gene_fun_units <- data.frame(coding = NA, fiveUTR = NA, interngenic = NA, intron = NA, promoter = NA, threeUTR = NA)

gene_fun_units[1,] <- as.numeric(table(vcf$genomic_func))



for (i in 1:length(sample_ids)){
  
  print(i)
  
  sample_id <- sample_ids[i]
  
  vcf <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF/annotated-variants/lung/",
                                sample_id, "/annotated-", sample_id, ".txt"),
                  header = T, sep = "\t", stringsAsFactors = F)
  
  gene_fun_units[i,] <- as.numeric(table(vcf$genomic_func))
  
}

# saveRDS(gene_fun_units, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/genic-units/func_ann_lung.rds")


png( filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/genic-units/mutation-distribution.png")
boxplot(log10(gene_fun_units), main = "Mutation occurrence in different genomic units", ylab = "log10 number of somatic mutations", xlab = "Genic functional units")
dev.off()
