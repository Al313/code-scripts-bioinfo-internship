#In this script I will try to reproduce fig. 3b of Akdemir paper. More specifically I will use the lung cohort and plot the proportion of SBS4 in active chrmatin regions to inactive regions.

library(tidyverse)
library(magrittr)
library(ggplot2)


options(scipen=999)

chromatin_bin_bed <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/data/reproducing-fig-3a/making-bed-files/chromatin-status-bins.bed.gz",
                              header = T, stringsAsFactors = F, sep = "\t")
chromatin_bin_bed %<>% mutate(span = chromatin_bin_bed$end - chromatin_bin_bed$start)
active_inactive_genome <- sum(chromatin_bin_bed$span[chromatin_bin_bed$chromatin_stat == "Active"])/sum(chromatin_bin_bed$span[chromatin_bin_bed$chromatin_stat == "Inactive"])

# Select the cohort that you want to use
sample_ids <- scan("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/sample-ids/annotated-breast-ids.txt", what = character())


# Select the SBS of interest and change 11 occurrances (change the case for 3 of them)
myMat <- NULL
j <- 0
for (i in 1:length(sample_ids)){
  
  
  sample_id <- sample_ids[i]
  vcf <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF/annotated-variants/breast/",
                                sample_id, "/annotated-", sample_id, ".txt"),
                  header = T, sep = "\t", stringsAsFactors = F)
  
  
  mut_sig_sum <- as.data.frame(table(vcf$mut_sign))
  
  if ("SBS2" %in% mut_sig_sum$Var1){
    print(sample_id)
    j <- j + 1
    
    mut_load <- nrow(vcf)
    mut_sig_sum %<>% mutate(contribution = mut_sig_sum$Freq/sum(mut_sig_sum$Freq))
    
    log2_sbs2_contribution <- log2(mut_sig_sum$contribution[mut_sig_sum$Var1 == "SBS2"]*100)
    
    
    
    chromatin_status_bins_sum <- as.data.frame(table(vcf$chromatin_status_bins))
    chromatin_status_bins_sum %<>% mutate(proportion = chromatin_status_bins_sum$Freq/sum(chromatin_status_bins_sum$Freq))
    
    
    parsed_vcf <- vcf[vcf$mut_sign == "SBS2" & vcf$chromatin_status_bins %in% c("Active", "Inactive"),]
    
    parsed_vcf_a <- as.numeric(parsed_vcf$mut_sign_score)
    mean_mut_score <- mean(parsed_vcf_a, na.rm = TRUE)
    
    
    active_inactive_sbs2 <- nrow(parsed_vcf[parsed_vcf$chromatin_status_bins == "Active",])/nrow(parsed_vcf[parsed_vcf$chromatin_status_bins == "Inactive",])
    
    active_inactive_sbs2_normalized <- active_inactive_sbs2/active_inactive_genome
    
    if (is.null(myMat)){
      print("null")
      myMat <- matrix(c(log2_sbs2_contribution, active_inactive_sbs2_normalized, mean_mut_score, mut_load), nrow = 1, dimnames = list(c(sample_id), c("log2_contribution", "chrom_ratio", "average_mut_score", "snv_load")))
      
    } else{
      myMat %<>% rbind(c(log2_sbs2_contribution ,active_inactive_sbs2_normalized, mean_mut_score,mut_load))
      rownames(myMat)[j] <- sample_ids[i]
    }
    
  }
  
}



myDf <- as.data.frame(myMat)
head(myDf)


myDf <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/chromatin-status/breast_mut_chromatin.rds")
saveRDS(myDf, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/chromatin-status/breast_mut_chromatin.rds")

#png(filename="/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/chromatin-status/sbs3(3).png")

ggplot(myDf, aes(x=chrom_ratio, y=log2_contribution)) +
  geom_point(aes(colour=snv_load)) +
  scale_color_gradient(low="blue", high="red", na.value = "grey50") +
  theme_classic() + 
  geom_smooth(method = lm, se = F) +
  labs(y="log2 (contribution to mutation load) (%)", 
       x="Normalized Mutation Ratio \n (active/inactive domain)", 
       title="Distribution of SBS4 mutations \n based on chromatin status \n \n \n", 
       caption = "HMF Lung Cohort (562 samples)") +
  theme(plot.title = element_text(hjust=0.5, face = "bold", size = 15),
        axis.text = element_text(face = "italic", size = 10, color = "black"),
        axis.title = element_text(face = "italic", size = 12, color = "black"))

#dev.off()
