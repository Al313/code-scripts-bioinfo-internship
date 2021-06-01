# Let's explore the signature analysis of data a bit in regards to their ber/ner status!

library("dplyr")
library("ggplot2")
library("tidyr")


if (dir.exists("/hpc/cuppen/")){
  metadata <- read.csv(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/labelled_metadata2.tsv", sep = "\t", header = T, stringsAsFactors = F)
  snv_count_dataset <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/metadata_with_snv_counts.rds")
  sample_dir <- "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF/annotated-variants/cancer-type/"
} else {
  metadata <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/labelled_metadata2.tsv", sep = "\t", header = T, stringsAsFactors = F)
  snv_count_dataset <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/metadata_with_snv_counts.rds")
  sample_dir <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF/annotated-variants/cancer-type/"
}


# getting the snv_counts in the main data matrix
metadata <- left_join(metadata, snv_count_dataset[,c("sampleId", "snv_count")], by = "sampleId")


# ## First let's take all the mutational signatures that there are
# 
# mut_sig <- vector()
# for (i in 1:nrow(metadata)){
#   message(i)
#   vcf <- read.csv(file = paste0(sample_dir, metadata$primaryTumorLocation[i], "/", metadata$sampleId[i], "/annotated-", metadata$sampleId[i], ".txt"), sep = "\t", header = T, stringsAsFactors = F)
#   container <- names(table(vcf$mut_sign))
#   mut_sig <- union(mut_sig, container)
# 
# }
# 
# 
# saveRDS(mut_sig, file = "./mut_sigs.rds")

if (dir.exists("/hpc/cuppen/")){
  mut_sigs <- readRDS(file = "/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-likelihood-score/mut_sigs.rds")
} else {
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


# # getting the mut sig counts in one matrix
# signature_matrix <- data.frame(sampleId = NA, primaryTumorLoc = NA, pathway = NA, sub_pathway = NA, hgnc_symbol = NA, snv_count = NA)
# signature_matrix[mut_sigs] <- NA
# 
# for (i in 1:nrow(metadata)){
#   message(i)
#   signature_matrix[i,c(1,2,3,4,5,6)] <- metadata[i,c(2,7,29,30,31,32)]
#   vcf <- read.csv(file = paste0(sample_dir, metadata$primaryTumorLocation[i], "/", metadata$sampleId[i], "/annotated-", metadata$sampleId[i], ".txt"), sep = "\t", header = T, stringsAsFactors = F)
# 
#   nn <- names(table(vcf$mut_sign))
#   vv <- as.vector(table(vcf$mut_sign))
# 
#   signature_matrix[i,nn] <- vv
# 
# }
# 
# 
# saveRDS(signature_matrix, file = "./mut_sigs_matrix_likelihood_score.rds")





signature_matrix_likelihood_approach <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-likelihood-score/mut_sigs_matrix_likelihood_score.rds")


signature_matrix_likelihood_approach[,7:60] <- signature_matrix_likelihood_approach[,7:60]/signature_matrix_likelihood_approach[,6]


# getting the mean contribution for each mut sig for each sample group _ for those samples that have the mutational signature
my_df_linkelihood_mut_sig <- data.frame(pathway = NA)
my_df_linkelihood_mut_sig[mut_sigs] <- NA

i <- 0

for (pathway in unique(signature_matrix_likelihood_approach$pathway)[c(1,2,4,5,6)]){
  i <- i + 1
  my_df_linkelihood_mut_sig[i,"pathway"] <- pathway

  for (mut_sig in mut_sigs){
    print(mut_sig)
    my_df_linkelihood_mut_sig[i,mut_sig] <- mean(signature_matrix_likelihood_approach[signature_matrix_likelihood_approach$pathway == pathway, mut_sig], na.rm = T)
    

  }

}


my_df_likelihood_mut_sig_tibble <- gather_(my_df_linkelihood_mut_sig,
                                           key_col = "mut_sigs",
                                           value_col = "value",
                                           gather_cols = colnames(my_df_linkelihood_mut_sig[,2:55]))

my_df_likelihood_mut_sig_tibble[,1] <- factor(my_df_likelihood_mut_sig_tibble[,1], levels = c("WT","BER","NER","POL","MMR"))
my_df_likelihood_mut_sig_tibble[,2] <- factor(my_df_likelihood_mut_sig_tibble[,2], levels = sig_order)
my_df_likelihood_mut_sig_tibble[,3] <- as.numeric(my_df_likelihood_mut_sig_tibble[,3])


my_df_likelihood_mut_sig_tibble$value <- my_df_likelihood_mut_sig_tibble$value * 100

# removing bars that would be present in less than 3 samples in the final plot!
# First we need to get the sample counts
signature_matrix <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-likelihood-score/mut_sigs_matrix_likelihood_score.rds")
sample_count <- data.frame(mut_sigs = NA, pathway = NA, nr_sample = NA)
i <- 0


for (sig in mut_sigs){
  for (pathway in c("WT", "BER", "MMR", "POL","NER")){
    i <- i + 1
    sample_count[i,"mut_sigs"] <- sig
    sample_count$pathway[i] <- pathway
    sample_count$nr_sample[i] <- sum(signature_matrix[,"pathway"] == pathway & !is.na(signature_matrix[,sig]))

  }
}


# This data frame contains the number of sample for each bar and I want to use it in the geom_text()
sample_count$value <- my_df_likelihood_mut_sig_tibble$value + 3
sample_count$mut_sigs <- factor(sample_count$mut_sigs, sig_order)





# demo data for plotting
my_df_likelihood_mut_sig_tibble1 <- my_df_likelihood_mut_sig_tibble[my_df_likelihood_mut_sig_tibble$mut_sigs %in% c("SBS1","SBS2","SBS10a","SBS15","SBS44"),]
sample_count1 <- sample_count[sample_count$mut_sigs %in% c("SBS1","SBS2","SBS10a","SBS15","SBS44"),]

# getting the barplots that have less than 3 samples
data_to_be_removed <- sample_count[sample_count$nr_sample > 0 & sample_count$nr_sample < 3, ]

i <- 1

for (i in 1:nrow(data_to_be_removed)) {
  my_df_likelihood_mut_sig_tibble$value[my_df_likelihood_mut_sig_tibble$mut_sigs == data_to_be_removed$mut_sigs[i] & my_df_likelihood_mut_sig_tibble$pathway == data_to_be_removed$pathway[i]] <- NA
}

i <- 1

for (i in 1:nrow(data_to_be_removed)) {
  sample_count$value[sample_count$mut_sigs == data_to_be_removed$mut_sigs[i] & sample_count$pathway == data_to_be_removed$pathway[i]] <- 0
}


w <- ggplot(my_df_likelihood_mut_sig_tibble, aes(x=pathway, y=value)) +
  geom_bar(stat = "identity", width = 0.7) + facet_wrap(~ mut_sigs, ncol = 5) +
  geom_text(data = sample_count, aes(x = pathway, y = value, label = nr_sample, fill = NULL),nudge_y=0.02, size = 3, color = 'red') +
  theme_bw() +
  ggtitle("Mean contribution of each mutational signature (likelihood-score approach)") +
  labs(y = "value %") +
  theme(plot.title = element_text(colour = "black",  face = "bold.italic", family = "Helvetica", size = 15, hjust = 0.5)) +
  theme(strip.text = element_text(size=10, face="bold", color="darkblue")) +
  theme(strip.background = element_rect(fill="lightblue", size=1, color="black"))

w

ggsave(filename = "./likelihood-mut-sig-contribution-corrected.pdf", plot = w, width = 50, height = 110, units = "cm")


if (dir.exists("/hpc/cuppen/")){
  ggsave(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-likelihood-score/likelihood-mut-sig-contribution-corrected.pdf",
         plot = w, width = 50, height = 110, units = "cm")
} else {
  ggsave(filename = "./likelihood-mut-sig-contribution-corrected.pdf", plot = w, width = 50, height = 110, units = "cm")
  ggsave(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-likelihood-score/likelihood-mut-sig-contribution-corrected.pdf", plot = w, width = 50, height = 110, units = "cm")
}





######## pole_linx




pole_linx <- scan(file = '/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/sample-ids/driver-genes/pole/pole-sample-ids.txt', what = character())

signature_matrix_likelihood_approach$pole_linx <- signature_matrix_likelihood_approach$sampleId %in% pole_linx

signature_matrix_likelihood_approach <- cbind(signature_matrix_likelihood_approach[,1:6], pole_linx = signature_matrix_likelihood_approach[,61], signature_matrix_likelihood_approach[,7:60])

signature_matrix_likelihood_approach[signature_matrix_likelihood_approach$pole_linx, "pole_linx"] <- "POLE"
signature_matrix_likelihood_approach[signature_matrix_likelihood_approach$pole_linx == "FALSE", "pole_linx"] <- "Non-POLE"


signature_matrix_likelihood_approach$sampleId[signature_matrix_likelihood_approach$pole_linx == "POLE" & !is.na(signature_matrix_likelihood_approach$SBS10a)]



my_df_likelihood_mut_sig_linx <- data.frame(pathway = NA)
my_df_likelihood_mut_sig_linx[mut_sigs] <- NA

i <- 0

for (pathway in unique(signature_matrix_likelihood_approach$pole_linx)){

  i <- i + 1
  my_df_likelihood_mut_sig_linx[i,"pathway"] <- pathway

  for (mut_sig in mut_sigs){
    print(mut_sig)
    my_df_likelihood_mut_sig_linx[i,mut_sig] <- mean(signature_matrix_likelihood_approach[signature_matrix_likelihood_approach$pole_linx == pathway, mut_sig], na.rm = T)

  }

}

my_df_likelihood_mut_sig_linx_tibble <- gather_(my_df_likelihood_mut_sig_linx,
                                                key_col = "mut_sigs",
                                                value_col = "value",
                                                gather_cols = colnames(my_df_likelihood_mut_sig_linx[,2:55]))

my_df_likelihood_mut_sig_linx_tibble[,1] <- factor(my_df_likelihood_mut_sig_linx_tibble[,1], levels = c("Non-POLE", "POLE"))
my_df_likelihood_mut_sig_linx_tibble[,2] <- factor(my_df_likelihood_mut_sig_linx_tibble[,2], levels = sig_order)
my_df_likelihood_mut_sig_linx_tibble[,3] <- as.numeric(my_df_likelihood_mut_sig_linx_tibble[,3])


my_df_likelihood_mut_sig_linx_tibble$value <- my_df_likelihood_mut_sig_linx_tibble$value * 100

# removing bars that would trpresent less than 3 samples in the final plot!
# First we need to get the sample counts
sample_count <- data.frame(mut_sigs = NA, pathway = NA, nr_sample = NA)
i <- 0


for (sig in mut_sigs){
  for (pathway in c("Non-POLE", "POLE")){
    i <- i + 1
    sample_count[i,"mut_sigs"] <- sig
    sample_count$pathway[i] <- pathway
    sample_count$nr_sample[i] <- sum(signature_matrix_likelihood_approach[,"pole_linx"] == pathway & !is.na(signature_matrix_likelihood_approach[,sig]))

  }
}


# This data frame contains the number of sample for each bar and I want to use it in the geom_text()
sample_count$value <- my_df_likelihood_mut_sig_linx_tibble$value + 1
sample_count$mut_sigs <- factor(sample_count$mut_sigs, sig_order)




# getting the barplots that have less than 3 samples
data_to_be_removed <- sample_count[sample_count$nr_sample > 0 & sample_count$nr_sample < 3, ]

i <- 1

for (i in 1:nrow(data_to_be_removed)) {
  my_df_likelihood_mut_sig_linx_tibble$value[my_df_likelihood_mut_sig_linx_tibble$mut_sigs == data_to_be_removed$mut_sigs[i] & my_df_likelihood_mut_sig_linx_tibble$pathway == data_to_be_removed$pathway[i]] <- NA
}

i <- 1

for (i in 1:nrow(data_to_be_removed)) {
  sample_count$value[sample_count$mut_sigs == data_to_be_removed$mut_sigs[i] & sample_count$pathway == data_to_be_removed$pathway[i]] <- 0
}



e <- ggplot(my_df_likelihood_mut_sig_linx_tibble, aes(x=pathway, y=value)) +
  geom_bar(stat = "identity", width = 0.7) + facet_wrap(~ mut_sigs, ncol = 5) +
  geom_text(data = sample_count, aes(x = pathway, y = value, label = nr_sample, fill = NULL),nudge_y=0.02, size = 3, color = 'red') +
  theme_bw() +
  ggtitle("Mean contribution of each mutational signature (likelihood-score approach) for pole_linx labeled samples") +
  labs(y = "Mean Contribution (%)") +
  theme(plot.title = element_text(colour = "black",  face = "bold.italic", family = "Helvetica", size = 15, hjust = 0.5)) +
  theme(strip.text = element_text(size=10, face="bold", color="darkblue")) +
  theme(strip.background = element_rect(fill="lightblue", size=1, color="black"))

e

if (dir.exists("/hpc/cuppen/")){
  ggsave(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-likelihood-score/likelihood-mut-sig-contribution_pole_linx.pdf",
         plot = e, width = 50, height = 110, units = "cm")
} else {
  ggsave(filename = "./likelihood-mut-sig-contribution_pole_linx.pdf", plot = e, width = 50, height = 110, units = "cm")
  ggsave(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-likelihood-score/likelihood-mut-sig-contribution_pole_linx.pdf", plot = e, width = 50, height = 110, units = "cm")
}




## Making the heatmap

library(pheatmap)
library(RColorBrewer)

nameList <- list(rownames = mut_sigs, colnames = c("BER/CON", "MMR/CON", "POL/CON", "NER/CON"))
mymatrix <- matrix(nrow = length(mut_sigs), ncol = 4, dimnames = nameList)


for (mut_sig in mut_sigs){
  for (i in 1:ncol(mymatrix)){
    pathway <- rownames(my_df_linkelihood_mut_sig)[i]
    j <- i+1
    mymatrix[mut_sig,i] <- log2(my_df_linkelihood_mut_sig[j,mut_sig]/my_df_linkelihood_mut_sig[1,mut_sig])
  }
}

mymatrix[is.na(mymatrix)] <- 0

myBreaks <- seq(from = -3.75, to = 3.75, by = 0.75)
myColors <- colorRampPalette(c("red", "white", "blue"))(length(myBreaks))

myColors[5] <- "#FFFFFF"



sample_count <- data.frame(pathway = NA)
sample_count[mut_sigs] <- NA

i <- 0

for (pathway in unique(signature_matrix_likelihood_approach$pathway)[c(1,2,4,5,6)]){
  i <- i + 1
  sample_count[i,"pathway"] <- pathway
  
  for (mut_sig in mut_sigs){
    print(mut_sig)
    sample_count[i,mut_sig] <- sum(!is.na(signature_matrix_likelihood_approach[signature_matrix_likelihood_approach$pathway == pathway, mut_sig]))
    
    
  }
  
}

sample_count_tibble <- gather_(sample_count, key_col = "mut_sigs", value_col = "nr_sample", gather_cols = colnames(sample_count[,2:55]))



# getting the barplots that have less than 3 samples
data_to_be_removed <- sample_count_tibble[sample_count_tibble$nr_sample > 0 & sample_count_tibble$nr_sample < 3, ]

i <- 1

for (i in 1:nrow(data_to_be_removed)) {
  
  if (data_to_be_removed$pathway[i] == "BER"){
    j <- "BER/CON"
  } else if (data_to_be_removed$pathway[i] == "NER"){
    j <- "NER/CON"
  } else if (data_to_be_removed$pathway[i] == "POL"){
    j <- "POL/CON"
  } else if (data_to_be_removed$pathway[i] == "MMR"){
    j <- "MMR/CON"
  }
  
  mymatrix[data_to_be_removed$mut_sigs[i],j] <- 0
}




png(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-likelihood-score/mut_sig_heat_map_corrected.png")
pheatmap(mymatrix, breaks = myBreaks, color = myColors, cluster_rows=FALSE, cluster_cols=FALSE, fontsize_row = 5, legend_breaks = myBreaks, legend_labels = c("", as.character(myBreaks)[2:10], "Log2\n"), main = "Log2 chnage of mutational signature contribution \n (likelihood-score method)")
dev.off()
while (!is.null(dev.list())) dev.off()



