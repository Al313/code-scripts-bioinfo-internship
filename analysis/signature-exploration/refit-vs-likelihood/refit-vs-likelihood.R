# ============| Correlation plot for contribution of SBS10a based on refit and linkelihood methods |===================



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


refit_sig_cont <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/annotation-beds/refit-signatures-contribution.txt.gz", sep = "\t", header = T, stringsAsFactors = F)
refit_sig_cont <- refit_sig_cont[,1:53]
refit_sig_cont <- refit_sig_cont[rownames(refit_sig_cont) %in% metadata$sampleId,]
refit_sig_cont <- cbind(sampleId = rownames(refit_sig_cont), refit_sig_cont)
rownames(refit_sig_cont) <- 1:nrow(refit_sig_cont)
# refit_sig_cont[refit_sig_cont$sampleId == "DRUP01010104T",]
refit_sig_cont <- refit_sig_cont[order(refit_sig_cont$sampleId),]



signature_matrix_likelihood_approach <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/sig-likelihood-score/mut_sigs_matrix_likelihood_score.rds")
# signature_matrix_likelihood_approach[signature_matrix_likelihood_approach$sampleId == "DRUP01010104T",]
signature_matrix_likelihood_approach <- signature_matrix_likelihood_approach[order(signature_matrix_likelihood_approach$sampleId),]


pole_linx <- scan(file = '/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/sample-ids/driver-genes/pole/pole-sample-ids.txt', what = character())

signature_matrix_likelihood_approach$pole_linx <- signature_matrix_likelihood_approach$sampleId %in% pole_linx
ncol(signature_matrix_likelihood_approach)
signature_matrix_likelihood_approach <- cbind(signature_matrix_likelihood_approach[,1:6], pole_linx = signature_matrix_likelihood_approach[,61], signature_matrix_likelihood_approach[,7:60])

signature_matrix_likelihood_approach[signature_matrix_likelihood_approach$pole_linx, "pole_linx"] <- "POLE"
signature_matrix_likelihood_approach[signature_matrix_likelihood_approach$pole_linx == "FALSE", "pole_linx"] <- "Non-POLE"








SBS10_active_samples_shared <- intersect(signature_matrix_likelihood_approach[!is.na(signature_matrix_likelihood_approach$SBS10a),"sampleId"], refit_sig_cont$sampleId[refit_sig_cont$SBS10a != 0])




refit_sbs10a_values <- refit_sig_cont[refit_sig_cont$sampleId %in% SBS10_active_samples_shared, "SBS10a"]
likelihood_sbs10a_values <- signature_matrix_likelihood_approach[signature_matrix_likelihood_approach$sampleId %in% SBS10_active_samples_shared, "SBS10a"]
pathway <- signature_matrix_likelihood_approach[signature_matrix_likelihood_approach$sampleId %in% SBS10_active_samples_shared, "pathway"]
linx_pathway <- signature_matrix_likelihood_approach[signature_matrix_likelihood_approach$sampleId %in% SBS10_active_samples_shared, "pole_linx"]

my_df_sbs10a <- data.frame(sampleId = SBS10_active_samples_shared, refit_values = refit_sbs10a_values, likelihood_values = likelihood_sbs10a_values, 
                           pathway = pathway, linx_pathway = linx_pathway)


minors <- seq(from = 0, to = 15, by = 2)
reg <- lm(formula = log(refit_sbs10a_values)~log(likelihood_sbs10a_values))


refit_vs_likelihood_plot_linx <- ggplot(my_df_sbs10a, aes(x = log(refit_sbs10a_values), y = log(likelihood_values), color = linx_pathway)) + 
  geom_point() + 
#  geom_text(label = my_df_sbs10a$sampleId) +
  theme_classic() +
  geom_hline(yintercept = minors, color = "grey80") +
  geom_vline(xintercept = minors, color = "grey80") +
  ggtitle(label = "Absolute contribution of SBS10a \n based on refit and likelihood approach") +
  theme(plot.title = element_text(colour = "black",  face = "bold.italic", family = "Helvetica", size = 15, hjust = 0.5)) +
  # geom_smooth(method='lm', formula= y~x, se = F) +
  scale_color_manual(values=c("#cc6699", "#0000ff"), name = "LINX Annotation")



if (dir.exists("/hpc/cuppen/")){
  ggsave(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/refit-vs-likelihood/refit-vslikelihood-linx.pdf", 
         plot = refit_vs_likelihood_plot_linx, width = 20, height = 20, units = "cm")
} else {
  ggsave(filename = "./refit-vslikelihood-linx.pdf", plot = refit_vs_likelihood_plot_linx, width = 20, height = 20, units = "cm")
  ggsave(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/refit-vs-likelihood/refit-vslikelihood-linx.pdf", plot = refit_vs_likelihood_plot_linx, width = 20, height = 20, units = "cm")
}


refit_vs_likelihood_plot_luan <- ggplot(my_df_sbs10a, aes(x = log(refit_sbs10a_values), y = log(likelihood_values), color = pathway)) + 
  geom_point(alpha = 0.65) + 
  theme_classic() +
  geom_hline(yintercept = minors, color = "grey80") +
  geom_vline(xintercept = minors, color = "grey80") +
  ggtitle(label = "Absolute contribution of SBS10a \n based on refit and likelihood approach") +
  theme(plot.title = element_text(colour = "black",  face = "bold.italic", family = "Helvetica", size = 15, hjust = 0.5)) +
  # geom_smooth(method='lm', formula= y~x, se = F) +
  scale_color_manual(values=c("#00ff00", "#0000ff", "#ffff00", "#cc6699", "#ff4da6"), name = "Luan Annotation")





if (dir.exists("/hpc/cuppen/")){
  ggsave(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/refit-vs-likelihood/refit-vslikelihood-luan.pdf", 
         plot = refit_vs_likelihood_plot_luan, width = 20, height = 20, units = "cm")
} else {
  ggsave(filename = "./refit-vslikelihood-luan.pdf", plot = refit_vs_likelihood_plot_luan, width = 20, height = 20, units = "cm")
  ggsave(filename = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/signature-exploration/refit-vs-likelihood/refit-vslikelihood-luan.pdf", plot = refit_vs_likelihood_plot_luan, width = 20, height = 20, units = "cm")
}


































