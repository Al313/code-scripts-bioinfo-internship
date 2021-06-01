# Loading libraries
library(magrittr)
library(dplyr)

# Select the cohort that you want to use
sample_ids <- scan("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/sample-ids/annotated-lung-ids.txt", what = character())

# sample_id <- sample_ids[1]


tobacco_df <- data.frame(c2a_unknown = NA, c2a_trp = NA, c2a_untrp = NA, g2t_unknown = NA, g2t_trp = NA, g2t_untrp = NA, c2t_unknown = NA, c2t_trp = NA, c2t_untrp = NA, g2a_unknown = NA, g2a_trp = NA, g2a_untrp = NA, tot_unknown = NA, tot_trp = NA, tot_untrp = NA)
non_tobacco_df <- data.frame(c2a_unknown = NA, c2a_trp = NA, c2a_untrp = NA, g2t_unknown = NA, g2t_trp = NA, g2t_untrp = NA, c2t_unknown = NA, c2t_trp = NA, c2t_untrp = NA, g2a_unknown = NA, g2a_trp = NA, g2a_untrp = NA, tot_unknown = NA, tot_trp = NA, tot_untrp = NA)

for (i in 1:length(sample_ids)){
  
  print(i)
  
  sample_id <- sample_ids[i]
  
  vcf <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF/annotated-variants/lung/",
                                sample_id, "/annotated-", sample_id, ".txt"),
                  header = T, sep = "\t", stringsAsFactors = F)
  
  mut_sig_sum <- as.data.frame(table(vcf$mut_sign))
  
  vcf_c2a <- vcf[vcf$REF == "C" & vcf$ALT == "A",]
  vcf_g2t <- vcf[vcf$REF == "G" & vcf$ALT == "T",]
  vcf_c2t <- vcf[vcf$REF == "C" & vcf$ALT == "T",]
  vcf_g2a <- vcf[vcf$REF == "G" & vcf$ALT == "A",]
  
  
  if ("SBS4" %in% mut_sig_sum$Var1){
    print(paste0(sample_id, " is smoking-associated!"))
    
    mut_sig_sum %<>% mutate(contribution = mut_sig_sum$Freq/sum(mut_sig_sum$Freq))
    
    if (mut_sig_sum$contribution[mut_sig_sum$Var1 == "SBS4"] > 0.2) {
      
      
      c2a <- as.numeric(table(vcf_c2a$trp_str_ann))
      tobacco_df[i,1:3] <- c2a
      g2t <- as.numeric(table(vcf_g2t$trp_str_ann))
      tobacco_df[i,4:6] <- g2t
      c2t <- as.numeric(table(vcf_c2t$trp_str_ann))
      tobacco_df[i,7:9] <- c2t
      g2a <- as.numeric(table(vcf_g2a$trp_str_ann))
      tobacco_df[i,10:12] <- g2a
      
      tot <- as.numeric(table(vcf$trp_str_ann))
      tobacco_df[i,13:15] <- tot
      
      
    } 
  
  } else {
    
    print(paste0(sample_id, " is not smoking-associated!"))
          
    c2a <- as.numeric(table(vcf_c2a$trp_str_ann))
    non_tobacco_df[i,1:3] <- c2a
    g2t <- as.numeric(table(vcf_g2t$trp_str_ann))
    non_tobacco_df[i,4:6] <- g2t
    c2t <- as.numeric(table(vcf_c2t$trp_str_ann))
    non_tobacco_df[i,7:9] <- c2t
    g2a <- as.numeric(table(vcf_g2a$trp_str_ann))
    non_tobacco_df[i,10:12] <- g2a
    
    tot <- as.numeric(table(vcf$trp_str_ann))
    non_tobacco_df[i,13:15] <- tot
    
  }
  
}

head(tobacco_df)
head(non_tobacco_df)

tobacco_df <- tobacco_df[,1:15]


# saveRDS(tobacco_df, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/trp-ori/tobacco_df.rds")
# saveRDS(o=non_tobacco_df, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/trp-ori/non_tobacco_df.rds")

png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/trp-ori/trp_str_bias_pos_con_smoking.png")

# C TO A in smoking-associated cases

tobacco_df <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/trp-ori/tobacco_df.rds")

tobacco_df %<>% mutate(sum_c2a = c2a_trp + c2a_untrp, sum_tot = tot_trp + tot_untrp)
tobacco_df <- tobacco_df[!is.na(tobacco_df$c2a_unknown),]

tobacco_df_normalized_c2a <- tobacco_df[,c(2,3)]/tobacco_df[,16]
tobacco_df_normalized_c2a[,c(3,4)] <- tobacco_df[,c(14,15)]/tobacco_df[,17]

# tobacco_df_normalized <- tobacco_df[,c(2,3,14,15)]

tobacco_df_normalized_c2a[,1] <- tobacco_df_normalized_c2a[,1]/tobacco_df_normalized_c2a[,3]
tobacco_df_normalized_c2a[,2] <- tobacco_df_normalized_c2a[,2]/tobacco_df_normalized_c2a[,4]


# boxplot(tobacco_df_normalized_c2a[,c(1,2)])


# G TO T in smoking-associated cases

tobacco_df <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/trp-ori/tobacco_df.rds")

tobacco_df %<>% mutate(sum_g2t = g2t_trp + g2t_untrp, sum_tot = tot_trp + tot_untrp)
tobacco_df <- tobacco_df[!is.na(tobacco_df$g2t_unknown),]

# tobacco_df_normalized <- tobacco_df[,c(5,6,14,15)]

tobacco_df_normalized_g2t <- tobacco_df[,c(5,6)]/tobacco_df[,16]
tobacco_df_normalized_g2t[,c(3,4)] <- tobacco_df[,c(14,15)]/tobacco_df[,17]

tobacco_df_normalized_g2t[,1] <- tobacco_df_normalized_g2t[,1]/tobacco_df_normalized_g2t[,3]
tobacco_df_normalized_g2t[,2] <- tobacco_df_normalized_g2t[,2]/tobacco_df_normalized_g2t[,4]


# boxplot(tobacco_df_normalized_g2t[,c(1,2)])


# combined smoking
tobacco_df_normalized <- cbind(tobacco_df_normalized_c2a,tobacco_df_normalized_g2t)
boxplot(tobacco_df_normalized[,c(1,2,5,6)], main = "Transcription Strand Bias in Smoking-Assocaited Lung Cancers", ylab = "Observed/Expected", xlab = "Mutation Types")

dev.off()



png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/trp-ori/trp_str_bias_pos_con_non_smoking.png")


# C TO A in non-smoking-associated cases

non_tobacco_df <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/trp-ori/non_tobacco_df.rds")

non_tobacco_df %<>% mutate(sum_c2a = c2a_trp + c2a_untrp, sum_tot = tot_trp + tot_untrp)
non_tobacco_df <- non_tobacco_df[!is.na(non_tobacco_df$c2a_unknown),]

non_tobacco_df_normalized_c2a <- non_tobacco_df[,c(2,3)]/non_tobacco_df[,16]
non_tobacco_df_normalized_c2a[,c(3,4)] <- non_tobacco_df[,c(14,15)]/non_tobacco_df[,17]

non_tobacco_df_normalized_c2a[,1] <- non_tobacco_df_normalized_c2a[,1]/non_tobacco_df_normalized_c2a[,3]
non_tobacco_df_normalized_c2a[,2] <- non_tobacco_df_normalized_c2a[,2]/non_tobacco_df_normalized_c2a[,4]


# boxplot(non_tobacco_df_normalized_c2a[,c(1,2)])


# G TO T in non-smoking-associated cases

non_tobacco_df <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/trp-ori/non_tobacco_df.rds")

non_tobacco_df %<>% mutate(sum_g2t = g2t_trp + g2t_untrp, sum_tot = tot_trp + tot_untrp)
non_tobacco_df <- non_tobacco_df[!is.na(non_tobacco_df$g2t_unknown),]

non_tobacco_df_normalized_g2t <- non_tobacco_df[,c(5,6)]/non_tobacco_df[,16]
non_tobacco_df_normalized_g2t[,c(3,4)] <- non_tobacco_df[,c(14,15)]/non_tobacco_df[,17]

non_tobacco_df_normalized_g2t[,1] <- non_tobacco_df_normalized_g2t[,1]/non_tobacco_df_normalized_g2t[,3]
non_tobacco_df_normalized_g2t[,2] <- non_tobacco_df_normalized_g2t[,2]/non_tobacco_df_normalized_g2t[,4]


# boxplot(non_tobacco_df_normalized_g2t[,c(1,2)])


# combined non-smoking
non_tobacco_df_normalized <- cbind(non_tobacco_df_normalized_c2a,non_tobacco_df_normalized_g2t)
boxplot(non_tobacco_df_normalized[,c(1,2,5,6)], main = "Transcription Strand Bias in Not Smoking-Assocaited Lung Cancers", ylab = "Observed/Expected", xlab = "Mutation Types")

dev.off()



png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/trp-ori/trp_str_bias_neg_con_smoking.png")


# C TO T in smoking-associated cases

tobacco_df <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/trp-ori/tobacco_df.rds")

tobacco_df %<>% mutate(sum_c2t = c2t_trp + c2t_untrp, sum_tot = tot_trp + tot_untrp)
tobacco_df <- tobacco_df[!is.na(tobacco_df$c2t_unknown),]

tobacco_df_normalized_c2t <- tobacco_df[,c(8,9)]/tobacco_df[,16]
tobacco_df_normalized_c2t[,c(3,4)] <- tobacco_df[,c(14,15)]/tobacco_df[,17]

tobacco_df_normalized_c2t[,1] <- tobacco_df_normalized_c2t[,1]/tobacco_df_normalized_c2t[,3]
tobacco_df_normalized_c2t[,2] <- tobacco_df_normalized_c2t[,2]/tobacco_df_normalized_c2t[,4]


# boxplot(tobacco_df_normalized_c2t[,c(1,2)])



# G TO A in smoking-associated cases

tobacco_df <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/trp-ori/tobacco_df.rds")

tobacco_df %<>% mutate(sum_g2a = g2a_trp + g2a_untrp, sum_tot = tot_trp + tot_untrp)
tobacco_df <- tobacco_df[!is.na(tobacco_df$g2a_unknown),]

tobacco_df_normalized_g2a <- tobacco_df[,c(11,12)]/tobacco_df[,16]
tobacco_df_normalized_g2a[,c(3,4)] <- tobacco_df[,c(14,15)]/tobacco_df[,17]

tobacco_df_normalized_g2a[,1] <- tobacco_df_normalized_g2a[,1]/tobacco_df_normalized_g2a[,3]
tobacco_df_normalized_g2a[,2] <- tobacco_df_normalized_g2a[,2]/tobacco_df_normalized_g2a[,4]


# boxplot(tobacco_df_normalized_g2a[,c(1,2)])


# Combined smoking
tobacco_df_normalized <- cbind(tobacco_df_normalized_c2t,tobacco_df_normalized_g2a)
boxplot(tobacco_df_normalized[,c(1,2,5,6)], main = "Transcription Strand Bias in Smoking-Assocaited Lung Cancers", ylab = "Observed/Expected", xlab = "Mutation Types")

dev.off()


png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/trp-ori/trp_str_bias_neg_con_non_smoking.png")


# C TO T in non-smoking-associated cases

non_tobacco_df <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/trp-ori/non_tobacco_df.rds")

non_tobacco_df %<>% mutate(sum_c2t = c2t_trp + c2t_untrp, sum_tot = tot_trp + tot_untrp)
non_tobacco_df <- non_tobacco_df[!is.na(non_tobacco_df$c2t_unknown),]

non_tobacco_df_normalized_c2t <- non_tobacco_df[,c(8,9)]/non_tobacco_df[,16]
non_tobacco_df_normalized_c2t[,c(3,4)] <- non_tobacco_df[,c(14,15)]/non_tobacco_df[,17]

non_tobacco_df_normalized_c2t[,1] <- non_tobacco_df_normalized_c2t[,1]/non_tobacco_df_normalized_c2t[,3]
non_tobacco_df_normalized_c2t[,2] <- non_tobacco_df_normalized_c2t[,2]/non_tobacco_df_normalized_c2t[,4]


# boxplot(non_tobacco_df_normalized_c2t[,c(1,2)])



# G TO A in non-smoking-associated cases

non_tobacco_df <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/trp-ori/non_tobacco_df.rds")

non_tobacco_df %<>% mutate(sum_g2a = g2a_trp + g2a_untrp, sum_tot = tot_trp + tot_untrp)
non_tobacco_df <- non_tobacco_df[!is.na(non_tobacco_df$g2a_unknown),]

non_tobacco_df_normalized_g2a <- non_tobacco_df[,c(11,12)]/non_tobacco_df[,16]
non_tobacco_df_normalized_g2a[,c(3,4)] <- non_tobacco_df[,c(14,15)]/non_tobacco_df[,17]

non_tobacco_df_normalized_g2a[,1] <- non_tobacco_df_normalized_g2a[,1]/non_tobacco_df_normalized_g2a[,3]
non_tobacco_df_normalized_g2a[,2] <- non_tobacco_df_normalized_g2a[,2]/non_tobacco_df_normalized_g2a[,4]


# boxplot(non_tobacco_df_normalized_g2a[,c(1,2)])


# combined non-smoking
non_tobacco_df_normalized <- cbind(non_tobacco_df_normalized_c2t,non_tobacco_df_normalized_g2a)
boxplot(non_tobacco_df_normalized[,c(1,2,5,6)], main = "Transcription Strand Bias in Not Smoking-Assocaited Lung Cancers", ylab = "Observed/Expected", xlab = "Mutation Types")

dev.off()




#

sample_id <- sample_ids[1]


vcf <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF/annotated-variants/lung/",
                              sample_id, "/annotated-", sample_id, ".txt"),
                header = T, sep = "\t", stringsAsFactors = F)

head(vcf[vcf$REF == "G" | vcf$REF == "A",], n = 100)
vcf <- vcf[vcf$REF == "G" | vcf$REF == "A",] 

vcf[i,]$REF == "G" | vcf[i,]$REF == "A"
i <- 1
for (i in 1:nrow(vcf)){
  
  if (vcf[i,]$REF == "G" | vcf[i,]$REF == "A"){
    if (vcf[i,]$trp_str_ann == "transcribed"){
      vcf[i,]$trp_str_ann <- "untranscribed"
    } else if (vcf[i,]$trp_str_ann == "untranscribed"){
      vcf[i,]$trp_str_ann <- "transcribed"
    }
    
  }
  
}

index_for_trp <- which(vcf$REF == "G" | vcf$REF == "A")
head(vcf[index_for_transcribed,])
index_for_transcribed <- which((vcf$REF == "A" | vcf$REF == "G") & vcf$trp_str_ann == "transcribed")
index_for_untranscribed <- which((vcf$REF == "A" | vcf$REF == "G") & vcf$trp_str_ann == "untranscribed")


vcf[index_for_transcribed,]$trp_str_ann <- "untranscribed"
vcf[index_for_untranscribed,]$trp_str_ann <- "transcribed"


mut_sig_sum <- as.data.frame(table(vcf$mut_sign))
vcf_c2a <- vcf[vcf$REF == "C" & vcf$ALT == "A",]
vcf_g2t <- vcf[vcf$REF == "G" & vcf$ALT == "T",]
mut_sig_sum %<>% mutate(contribution = mut_sig_sum$Freq/sum(mut_sig_sum$Freq))
table(vcf_g2t$trp_str_ann)[2]/table(vcf_g2t$trp_str_ann)[3]
table(vcf_c2a$trp_str_ann)[2]/table(vcf_c2a$trp_str_ann)[3]
table(vcf$trp_str_ann)[2]/table(vcf$trp_str_ann)[3]
table(vcf_g2t$trp_str_ann)
table(vcf_c2a$trp_str_ann)
table(vcf$trp_str_ann)


tobacco_df <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/trp-ori/tobacco_df.rds")

tobacco_df %<>% mutate(sum_g2t = g2t_trp + g2t_untrp, sum_tot = tot_trp + tot_untrp)
tobacco_df <- tobacco_df[!is.na(tobacco_df$g2t_unknown),]

tobacco_df_normalized <- tobacco_df[,c(5,6)]/tobacco_df[,16]
tobacco_df_normalized[,c(3,4)] <- tobacco_df[,c(14,15)]/tobacco_df[,17]

tobacco_df_normalized[,1] <- tobacco_df_normalized[,1]/tobacco_df_normalized[,3]
tobacco_df_normalized[,2] <- tobacco_df_normalized[,2]/tobacco_df_normalized[,4]


boxplot(tobacco_df_normalized[,c(1,2)])