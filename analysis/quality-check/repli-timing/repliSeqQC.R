library(tidyverse)
library(magrittr)

par(mfrow = c(1,2))
sss <- sss[-"DRUP01330020T"]


sample_ids <- scan("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/sample-ids/annotated-breast-ids.txt", what = character())
sample_ids <- sample_ids[sample_ids != "DRUP01330020T"]

# My data

repliBedMyData <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/data/warm-up/archive/repliseq/bed-union/repliseq_encode_mean_binned.bed.gz",
                      header = T, sep = "\t", stringsAsFactors = F)
repliSumMyData <- as.data.frame(table(repliBedMyData$bin))
repliSumMyData %<>%
  mutate(portion = repliSumMyData$Freq/sum(repliSumMyData$Freq))
repliSumMyData <- repliSumMyData[c(1,3,2),]


# Boxtel data

repliBedBoxtel <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/data/warm-up/archive/repliseq/boxtel-paper/boxtel_paper_all_RepliSeq_median.bed.gz",
                           header = T, sep = "\t", stringsAsFactors = F)
repliSumBoxtel <- as.data.frame(table(repliBedBoxtel$bin))
repliSumBoxtel %<>%
  mutate(portion = repliSumBoxtel$Freq/sum(repliSumBoxtel$Freq))
repliSumBoxtel <- repliSumBoxtel[c(1,3,2),]




repliTimingSum <- NULL
repliTimingSum2 <- NULL
for (i in 1:length(sample_ids)){
  print(i)
  sample_id <- sample_ids[i]
  vcf <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF/annotated-variants/breast",
                                        sample_id, "/annotated-", sample_id, ".txt"),
                          header = T, sep = "\t", stringsAsFactors = F)[,c(15,16)]
  
  vcfSum <- as.data.frame(table(vcf$rep_timing))
  vcfSum <- vcfSum[c(1,3,2),]
  vcfSum %<>%
    mutate(portion = vcfSum$Freq/sum(vcfSum$Freq))
  if (i == 1){
    repliTimingSum <- as.matrix(log2(vcfSum$portion/repliSumMyData$portion))
    rownames(repliTimingSum) <- vcfSum$Var1  
    }
  else {
  repliTimingSum <- cbind(repliTimingSum, as.matrix(log2(vcfSum$portion/repliSumMyData$portion)))
  
  }
  
  colnames(repliTimingSum)[i] <- sample_id
  
  
  
  
  vcfSum2 <- as.data.frame(table(vcf$rep_timing_boxtel))
  vcfSum2 <- vcfSum2[c(1,3,2),]
  vcfSum2 %<>%
    mutate(portion = vcfSum2$Freq/sum(vcfSum2$Freq))
  if (i == 1){
    repliTimingSum2 <- as.matrix(log2(vcfSum2$portion/repliSumBoxtel$portion))
    rownames(repliTimingSum2) <- vcfSum2$Var1  
  }
  else {
    repliTimingSum2 <- cbind(repliTimingSum2, as.matrix(log2(vcfSum2$portion/repliSumBoxtel$portion)))
    
  }
  
  colnames(repliTimingSum2)[i] <- sample_id
  
}


sdd <- apply(repliTimingSum, 1, sd, na.rm = T)/sqrt(ncol(repliTimingSum))
names(sdd) <- c("early", "mid", "late")

summ <- apply(repliTimingSum, 1, mean, na.rm = T)/1



mid <- barplot(summ, ylim = c(-0.6,0.6), main = "Replication Timing Quality Check \n (Breast cohort) \n My data", ylab = "Log2(observed/expected)", xlab = "bins", col = c("#144343", "#23cbcb", "#04fcfc"))

arrows(x0=mid, y0=summ-sdd, x1=mid, y1=summ+sdd, code=3, col="red", lwd=2, angle=90, length=0.05)

legend("topleft", legend = vcf1Sum$Var1[c(1,3,2)], fill = c("#144343", "#23cbcb", "#04fcfc"))




# Boxtel data


sdd2 <- apply(repliTimingSum2, 1, sd, na.rm = T)/sqrt(ncol(repliTimingSum2))
names(sdd2) <- c("early", "mid", "late")

summ2 <- apply(repliTimingSum2, 1, mean, na.rm = T)/1



mid2 <- barplot(summ2, ylim = c(-0.5,0.5), main = "Replication Timing Quality Check \n (Breast cohort) \n Boxtel paper", ylab = "Log2(observed/expected)", xlab = "bins", col = c("#144343", "#23cbcb", "#04fcfc"))

arrows(x0=mid2, y0=summ2-sdd2, x1=mid2, y1=summ2+sdd2, code=3, col="red", lwd=2, angle=90, length=0.05)

legend("topleft", legend = vcfSum2$Var1[c(1,3,2)], fill = c("#144343", "#23cbcb", "#04fcfc"))


