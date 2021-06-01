# ===========| trp-ori and mut-sig visualization |===============


# Loading packages
library(tidyr)
library(dplyr)
library(ggplot2)


# Reading in the data
if (dir.exists("/hpc/cuppen/")){
  metadata <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/labelled_metadata2.tsv", sep = "\t", header = T, stringsAsFactors = F)
  snv_count_dataset <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/metadata_with_snv_counts.rds")
  sample_dir <- "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF/annotated-variants/cancer-type/"
  mut_sigs <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-likelihood-score/mut_sigs.rds")
} else {
  metadata <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/labelled_metadata2.tsv", sep = "\t", header = T, stringsAsFactors = F)
  snv_count_dataset <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/metadata_with_snv_counts.rds")
  sample_dir <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF/annotated-variants/cancer-type/"
  mut_sigs <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-likelihood-score/mut_sigs.rds")
}


# having the mutation signatures in order
sig_order <- paste0("SBS", c("1", "2", "3", "4", "5", "6", "7a",
                             "7b", "7c", "7d", "8", "9", "10a", "10b",
                             "10c", "10d", "11", "12", "13", "14", "15", "16", "17a",
                             "17b", "18", "19", "20", "21", "22", "23",
                             "24", "25", "26", "28", "29", "30",
                             "31", "32", "33", "34", "35", "36", "37",
                             "38", "39", "40", "41", "42", "44", "84",
                             "85", "86", "87", "88", "89", "90", "91",
                             "92", "93", "94"))

mut_sigs <- intersect(sig_order, mut_sigs)

# getting the snv_counts in the main data matrix
metadata <- left_join(metadata, snv_count_dataset[,c("sampleId", "snv_count")], by = "sampleId")




# # getting the matrix
# 
# trp_ori <- NA
# signature_trp_ori_matrix <- data.frame(sampleId = NA, primaryTumorLoc = NA, pathway = NA, sub_pathway = NA, hgnc_symbol = NA, snv_count = NA, trp_ori = NA)
# signature_trp_ori_matrix[mut_sigs] <- NA
# k <- 1
# 
# for (i in 1:nrow(metadata)){
# 
#   message(i)
# 
#   vcf <- read.csv(file = paste0(sample_dir, metadata$primaryTumorLocation[i], "/", metadata$sampleId[i], "/annotated-", metadata$sampleId[i], ".txt"), sep = "\t", header = T, stringsAsFactors = F)
# 
#   if (length(trp_ori)==1){
#     trp_ori <- names(table(vcf$trp_str_ann))
#   }
# 
#   for (j in trp_ori){
#     signature_trp_ori_matrix[k,c(1,2,3,4,5,6)] <- metadata[i,c(2,7,29,30,31,32)]
#     signature_trp_ori_matrix[k,7] <- j
# 
# 
#     nn <- names(table(vcf$mut_sign[vcf$trp_str_ann == j]))
#     vv <- as.vector(table(vcf$mut_sign[vcf$trp_str_ann == j]))
# 
#     signature_trp_ori_matrix[k,nn] <- vv
# 
#     k <- k + 1
#   }
# 
# }
# 
# 
# 
# saveRDS(signature_trp_ori_matrix, file = "./sig_trp_ori_matrix.rds")




if (dir.exists("/hpc/cuppen/")){
    signature_trp_ori_matrix <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-and-trp-ori/sig_trp_ori_matrix.rds")
} else {
  signature_trp_ori_matrix <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-and-trp-ori/sig_trp_ori_matrix.rds")
}


# removing bars that would trpresent less than 3 samples in the final plot!
# First we need to get the sample counts
signature_matrix <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-likelihood-score/mut_sigs_matrix_likelihood_score.rds")
sample_count <- data.frame(mut_sigs = NA, pathway = NA, nr_sample = NA)
i <- 0

for (sig in mut_sigs){
  for (pathway in c("WT", "BER", "POL", "NER", "MMR")){
    i <- i + 1
    sample_count[i,"mut_sigs"] <- sig
    sample_count$pathway[i] <- pathway
    sample_count$nr_sample[i] <- sum(signature_matrix[,"pathway"] == pathway & !is.na(signature_matrix[,sig]))

  }
}


# This data frame contains the number of sample for each bar and I want to use it in the geom_text()
sample_count$value <- 1
sample_count$mut_sigs <- factor(sample_count$mut_sigs, sig_order)




# converting absolute mutation counts to percent
options(scipen=999)
for (sample_id in metadata$sampleId){
  message(rownames(metadata[metadata$sampleId == sample_id,]))
  for (mut_sig in mut_sigs){


    signature_trp_ori_matrix[signature_trp_ori_matrix$sampleId == sample_id, mut_sig] <- signature_trp_ori_matrix[signature_trp_ori_matrix$sampleId == sample_id, mut_sig]/sum(signature_trp_ori_matrix[signature_trp_ori_matrix$sampleId == sample_id, mut_sig])

  }

}



# getting the mean contribution for each mut sig for each sample group

my_df_trp_ori <- data.frame(pathway = NA, trp_ori = NA)
my_df_trp_ori[mut_sigs] <- NA

i <- 0

for (trp_ori in unique(signature_trp_ori_matrix$trp_ori)){
  for (pathway in unique(signature_trp_ori_matrix$pathway)[c(1,2,4,5,6)]){

    i <- i + 1
    my_df_trp_ori[i,"pathway"] <- pathway
    my_df_trp_ori[i,"trp_ori"] <- trp_ori

    for (mut_sig in mut_sigs){
      print(mut_sig)
      my_df_trp_ori[i,mut_sig] <- mean(signature_trp_ori_matrix[signature_trp_ori_matrix$pathway == pathway & signature_trp_ori_matrix$trp_ori == trp_ori, mut_sig], na.rm = T)

    }


  }


}



# making a tibble
sig_trp_ori_tibble <- gather_(my_df_trp_ori,
                              key_col = "mut_sigs",
                              value_col = "value",
                              gather_cols = colnames(my_df_trp_ori[,3:56]))



# getting the correct data types (note that with specifying the levels such we also determine their order in the ggplot function)
sig_trp_ori_tibble[,1] <- factor(sig_trp_ori_tibble[,1], levels = c("WT","BER","NER","POL", "MMR"))
sig_trp_ori_tibble[,2] <- factor(sig_trp_ori_tibble[,2], levels = c("-", "transcribed", "untranscribed"))
sig_trp_ori_tibble[,3] <- factor(sig_trp_ori_tibble[,3], levels = sig_order)
sig_trp_ori_tibble[,4] <- as.numeric(sig_trp_ori_tibble[,4])



# demo data for plotting
sig_trp_ori_tibble1 <- sig_trp_ori_tibble[sig_trp_ori_tibble$mut_sigs %in% c("SBS1","SBS2","SBS10a","SBS15","SBS44"),]
sample_count1 <- sample_count[sample_count$mut_sigs %in% c("SBS1","SBS2","SBS10a","SBS15","SBS44"),]



# getting the barplots that have less than 3 samples
data_to_be_removed <- sample_count[sample_count$nr_sample > 0 & sample_count$nr_sample < 3, ]

i <- 1

for (i in 1:nrow(data_to_be_removed)) {
  sig_trp_ori_tibble$value[sig_trp_ori_tibble$mut_sigs == data_to_be_removed$mut_sigs[i] & sig_trp_ori_tibble$pathway == data_to_be_removed$pathway[i]] <- NA
}


# plotting

g <- ggplot(sig_trp_ori_tibble, aes(x=pathway, y=value, fill=trp_ori)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) + facet_wrap(~ mut_sigs, ncol = 5) +
  scale_fill_brewer("Transcription Strand Orientation", palette = "Pastel2") +
  geom_text(aes(label = paste0("%", formatC(value*100, digits = 3))), position = position_stack(vjust = 0.5), size = 3) +
  geom_text(data = sample_count, aes(x = pathway, y = value, label = nr_sample, fill = NULL),nudge_y=0.02, size = 3, color = 'red') +
  theme_bw() +
  ggtitle("Percentage of somatic mutations caused by each mutational signature \n with regard to transcriptional strand orientation") +
  theme(plot.title = element_text(colour = "black",  face = "bold.italic", family = "Helvetica", size = 25, hjust = 0.5)) +
  theme(strip.text = element_text(size=10, face="bold", color="darkblue")) +
  theme(strip.background = element_rect(fill="lightblue", size=1, color="black"))
g



if (dir.exists("/hpc/cuppen/")){
  ggsave(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-and-trp-ori/pathway-mutsig-trpori-corrected.pdf",
         plot = g, width = 50, height = 110, units = "cm")
} else {
  ggsave(filename = "pathway-mutsig-trpori-corrected.pdf", plot = g, width = 50, height = 110, units = "cm")
  ggsave(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-and-trp-ori/pathway-mutsig-trpori-corrected.pdf", plot = g, width = 50, height = 110, units = "cm")
}




