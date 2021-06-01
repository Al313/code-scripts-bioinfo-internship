# Loading packages

library("tidyr")
library("ggplot2")
library("dplyr")

if (dir.exists("/hpc/cuppen/")){
  refit_sig_cont <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-contribution-refit/refit-sig-contribution.txt.gz",
                             sep = "\t", stringsAsFactors = F, header = T)
  metadata <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/labelled_metadata2.tsv", 
                       sep = "\t", header = T, stringsAsFactors = F)
  snv_count_dataset <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/metadata_with_snv_counts.rds")
  sample_dir <- "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF/annotated-variants/cancer-type/"
  mut_sigs <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-likelihood-score/mut_sigs.rds")
} else {
  refit_sig_cont <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-contribution-refit/refit-sig-contribution.txt.gz",
                             sep = "\t", stringsAsFactors = F, header = T)
  metadata <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/labelled_metadata2.tsv", 
                       sep = "\t", header = T, stringsAsFactors = F)
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

# # Removing mutational signatures for other mutation types other than SBS
# 
# refit_sig_cont <- refit_sig_cont[order(rownames(refit_sig_cont)),]
# refit_sig_cont <- refit_sig_cont[,1:53]
# 
# # Keeping the sample names that are also in the metadata matrix (one sample from every one patient)
# refit_sig_cont <- refit_sig_cont[rownames(refit_sig_cont) %in% metadata$sampleId,]
# 
# 
# 
# # getting the snv_counts in the main data matrix
# metadata <- left_join(metadata, snv_count_dataset[,c("sampleId", "snv_count")], by = "sampleId")
# metadata <- metadata[order(metadata$sampleId),]
# 
# 
# # contribution of each signature to the overall mutational catalogue of each sample
# 
# refit_sig_cont_perc <- apply(refit_sig_cont, MARGIN = 1, FUN = function(x1) x1/sum(x1))
# refit_sig_cont_perc <- as.data.frame(t(as.matrix(refit_sig_cont_perc)))
# 
# 
# # Getting the necessary metadata columns from the metadata df and including them in the mut-sig-cont-perc df
# refit_sig_cont_perc <- cbind(metadata[,c(2,7,29,30,31,32)], refit_sig_cont_perc)
# rownames(refit_sig_cont_perc) <- 1:nrow(refit_sig_cont_perc)
# 
# 
# saveRDS(refit_sig_cont_perc, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-contribution-refit/processed-refit-sig-contribution.rds")
# 





processed_refit_sig_cont <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-contribution-refit/processed-refit-sig-contribution.rds")





# for refit approach the SBS30 is not present
mut_sigs_for_refit <- mut_sigs[mut_sigs != "SBS33"]

# getting the mean contribution for each mut sig for each sample group
my_df_refit_mut_sig <- data.frame(pathway = NA)
my_df_refit_mut_sig[mut_sigs_for_refit] <- NA

i <- 0

for (pathway in unique(processed_refit_sig_cont$pathway)[c(1,2,4,5,6)]){

  i <- i + 1
  my_df_refit_mut_sig[i,"pathway"] <- pathway

  for (mut_sig in mut_sigs_for_refit){
    print(mut_sig)
    my_df_refit_mut_sig[i,mut_sig] <- mean(processed_refit_sig_cont[processed_refit_sig_cont$pathway == pathway, mut_sig], na.rm = T)

  }

}


my_df_refit_mut_sig_tibble <- gather_(my_df_refit_mut_sig,
                                      key_col = "mut_sigs",
                                      value_col = "value",
                                      gather_cols = colnames(my_df_refit_mut_sig[,2:54]))

my_df_refit_mut_sig_tibble[,1] <- factor(my_df_refit_mut_sig_tibble[,1], levels = c("WT","BER","NER","POL","MMR"))
my_df_refit_mut_sig_tibble[,2] <- factor(my_df_refit_mut_sig_tibble[,2], levels = mut_sigs_for_refit)
my_df_refit_mut_sig_tibble[,3] <- as.numeric(my_df_refit_mut_sig_tibble[,3])

my_df_refit_mut_sig_tibble$value <- my_df_refit_mut_sig_tibble$value * 100


# removing bars that would trpresent less than 3 samples in the final plot!
# First we need to get the sample counts
signature_matrix <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-likelihood-score/mut_sigs_matrix_likelihood_score.rds")
sample_count <- data.frame(mut_sigs = NA, pathway = NA, nr_sample = NA)
i <- 0



for (sig in mut_sigs_for_refit){
  for (pathway in c("WT", "BER", "MMR", "POL","NER")){
    i <- i + 1
    sample_count[i,"mut_sigs"] <- sig
    sample_count$pathway[i] <- pathway
    sample_count$nr_sample[i] <- sum(processed_refit_sig_cont[,"pathway"] == pathway & processed_refit_sig_cont[,sig] != 0)

  }
}


# This data frame contains the number of sample for each bar and I want to use it in the geom_text()
sample_count$value <- my_df_refit_mut_sig_tibble$value + 0.5
sample_count$mut_sigs <- factor(sample_count$mut_sigs, mut_sigs_for_refit)





# demo data for plotting
my_df_refit_mut_sig_tibble1 <- my_df_refit_mut_sig_tibble[my_df_refit_mut_sig_tibble$mut_sigs %in% c("SBS1","SBS2","SBS10a","SBS15","SBS44"),]
sample_count1 <- sample_count[sample_count$mut_sigs %in% c("SBS1","SBS2","SBS10a","SBS15","SBS44"),]



# getting the barplots that have less than 3 samples
data_to_be_removed <- sample_count[sample_count$nr_sample > 0 & sample_count$nr_sample < 3, ]

i <- 1

for (i in 1:nrow(data_to_be_removed)) {
  my_df_refit_mut_sig_tibble$value[my_df_refit_mut_sig_tibble$mut_sigs == data_to_be_removed$mut_sigs[i] & my_df_refit_mut_sig_tibble$pathway == data_to_be_removed$pathway[i]] <- NA
}

i <- 1

for (i in 1:nrow(data_to_be_removed)) {
  sample_count$value[sample_count$mut_sigs == data_to_be_removed$mut_sigs[i] & sample_count$pathway == data_to_be_removed$pathway[i]] <- 0
}



q <- ggplot(my_df_refit_mut_sig_tibble, aes(x=pathway, y=value)) +
  geom_bar(stat = "identity", width = 0.7) + facet_wrap(~ mut_sigs, ncol = 5) +
  geom_text(data = sample_count, aes(x = pathway, y = value, label = nr_sample, fill = NULL),nudge_y=0.02, size = 3, color = 'red') +
  theme_bw() +
  ggtitle("Mean contribution of each mutational signature (refit approach)") +
  labs(y = "value %") +
  theme(plot.title = element_text(colour = "black",  face = "bold.italic", family = "Helvetica", size = 15, hjust = 0.5)) +
  theme(strip.text = element_text(size=10, face="bold", color="darkblue")) +
  theme(strip.background = element_rect(fill="lightblue", size=1, color="black"))

q


if (dir.exists("/hpc/cuppen/")){
  ggsave(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-contribution-refit/refit-mut-sig-contribution-corrected.pdf",
         plot = q, width = 50, height = 110, units = "cm")
} else {
  ggsave(filename = "refit-mut-sig-contribution-corrected.pdf", plot = q, width = 50, height = 110, units = "cm")
  ggsave(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-contribution-refit/refit-mut-sig-contribution-corrected.pdf", plot = q, width = 50, height = 110, units = "cm")
}







######### pole_linx




pole_linx <- scan(file = '/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/sample-ids/driver-genes/pole/pole-sample-ids.txt', what = character())

processed_refit_sig_cont$pole_linx <- processed_refit_sig_cont$sampleId %in% pole_linx

processed_refit_sig_cont <- cbind(processed_refit_sig_cont[,1:6], pole_linx = processed_refit_sig_cont[,60], processed_refit_sig_cont[,7:59])

processed_refit_sig_cont[processed_refit_sig_cont$pole_linx, "pole_linx"] <- "POLE"
processed_refit_sig_cont[processed_refit_sig_cont$pole_linx == "FALSE", "pole_linx"] <- "Non-POLE"



# getting the mean contribution for each mut sig for each sample group
my_df_refit_mut_sig_linx <- data.frame(pathway = NA)
my_df_refit_mut_sig_linx[mut_sigs_for_refit] <- NA

i <- 0

for (pathway in unique(processed_refit_sig_cont$pole_linx)){

  i <- i + 1
  my_df_refit_mut_sig_linx[i,"pathway"] <- pathway

  for (mut_sig in mut_sigs_for_refit){
    print(mut_sig)
    my_df_refit_mut_sig_linx[i,mut_sig] <- mean(processed_refit_sig_cont[processed_refit_sig_cont$pole_linx == pathway, mut_sig], na.rm = T)

  }

}


my_df_refit_mut_sig_linx_tibble <- gather_(my_df_refit_mut_sig_linx,
                                      key_col = "mut_sigs",
                                      value_col = "value",
                                      gather_cols = colnames(my_df_refit_mut_sig_linx[,2:54]))

my_df_refit_mut_sig_linx_tibble[,1] <- factor(my_df_refit_mut_sig_linx_tibble[,1], levels = c("Non-POLE", "POLE"))
my_df_refit_mut_sig_linx_tibble[,2] <- factor(my_df_refit_mut_sig_linx_tibble[,2], levels = mut_sigs_for_refit)
my_df_refit_mut_sig_linx_tibble[,3] <- as.numeric(my_df_refit_mut_sig_linx_tibble[,3])

my_df_refit_mut_sig_linx_tibble$value <- my_df_refit_mut_sig_linx_tibble$value * 100



# removing bars that would trpresent less than 3 samples in the final plot!
# First we need to get the sample counts
sample_count <- data.frame(mut_sigs = NA, pathway = NA, nr_sample = NA)
i <- 0


for (sig in mut_sigs_for_refit){
  for (pathway in c("Non-POLE", "POLE")){
    i <- i + 1
    sample_count[i,"mut_sigs"] <- sig
    sample_count$pathway[i] <- pathway
    sample_count$nr_sample[i] <- sum(processed_refit_sig_cont[,"pole_linx"] == pathway & processed_refit_sig_cont[,sig] != 0)

  }
}


# This data frame contains the number of sample for each bar and I want to use it in the geom_text()
sample_count$value <- my_df_refit_mut_sig_linx_tibble$value + 1
sample_count$mut_sigs <- factor(sample_count$mut_sigs, mut_sigs_for_refit)




# getting the barplots that have less than 3 samples
data_to_be_removed <- sample_count[sample_count$nr_sample > 0 & sample_count$nr_sample < 3, ]

i <- 1

for (i in 1:nrow(data_to_be_removed)) {
  my_df_refit_mut_sig_linx_tibble$value[my_df_refit_mut_sig_linx_tibble$mut_sigs == data_to_be_removed$mut_sigs[i] & my_df_refit_mut_sig_linx_tibble$pathway == data_to_be_removed$pathway[i]] <- NA
}

i <- 1

for (i in 1:nrow(data_to_be_removed)) {
  sample_count$value[sample_count$mut_sigs == data_to_be_removed$mut_sigs[i] & sample_count$pathway == data_to_be_removed$pathway[i]] <- 0
}




r <- ggplot(my_df_refit_mut_sig_linx_tibble, aes(x=pathway, y=value)) +
  geom_bar(stat = "identity", width = 0.7) + facet_wrap(~ mut_sigs, ncol = 5) +
  geom_text(data = sample_count, aes(x = pathway, y = value, label = nr_sample, fill = NULL),nudge_y=0.02, size = 3, color = 'red') +
  theme_bw() +
  ggtitle("Mean contribution of each mutational signature (refit approach) for pole_linx labeled samples") +
  labs(y = "value %") +
  theme(plot.title = element_text(colour = "black",  face = "bold.italic", family = "Helvetica", size = 15, hjust = 0.5)) +
  theme(strip.text = element_text(size=10, face="bold", color="darkblue")) +
  theme(strip.background = element_rect(fill="lightblue", size=1, color="black"))



if (dir.exists("/hpc/cuppen/")){
  ggsave(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-contribution-refit/refit-mut-sig-contribution-pole-linx-corrected.pdf",
         plot = r, width = 50, height = 110, units = "cm")
} else {
  ggsave(filename = "refit-mut-sig-contribution-pole-linx-corrected.pdf", plot = r, width = 50, height = 110, units = "cm")
  ggsave(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-contribution-refit/refit-mut-sig-contribution-pole-linx-corrected.pdf", plot = r, width = 50, height = 110, units = "cm")
}

