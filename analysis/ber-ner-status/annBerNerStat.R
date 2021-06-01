# Loading packages ...
library(stringr)
library(magrittr)
library(dplyr)
library(RColorBrewer)

# Reading in the data
diplotypes <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/annotation-beds/diplotypes_germPonFilt.txt.gz",
         sep = "\t", header = T, stringsAsFactors = F)



# General inspection of the data

length(unique(diplotypes$sample))
sum(!str_detect(unique(diplotypes$sample), "^DO"))

## There are 4,784 HMF samples (all HMF samples available on HPC) and 1760 PCAWG samples (out of 1,774 samples available on HPC)



# get only the data for HMF samples

hmf_diplotypes <- diplotypes[str_detect(diplotypes$sample, "^DO", negate = T),]
rownames(hmf_diplotypes) <- 1:nrow(hmf_diplotypes)
length(unique(hmf_diplotypes$sample))



# Get a sense of data
str(hmf_diplotypes)
table(hmf_diplotypes[hmf_diplotypes$a1.max_score == 5 & hmf_diplotypes$a2.max_score == 5,]$diplotype_origin)
table(hmf_diplotypes$hgnc_symbol[hmf_diplotypes$diplotype_origin == "cnv_cnv"])
hmf_diplotypes[hmf_diplotypes$hgnc_symbol == "POLE",]

# Determine which genes are involved in BER and which genes in NER

specific_pathways <- c("BER", "GG-NER", "NER", "NER", "NER", "NER", "NER", "NER", "NER", "TC-NER", "TC-NER", "GLY-BER", "MMR", "GLY-BER", "MMR", "MMR", "GLY-BER", "GLY-BER", "GLY-BER", "GLY-BER", "GLY-BER", "MMR", "GLY-BER", "POL_PCNA", "MMR", "POL_B", "POL_D1", "POL_E", "POL_K", "GLY-BER", "GLY-BER", "GLY-BER", "NER", "GG-NER", "NER/BER")
collective_pathways <- c("BER", "NER", "NER", "NER", "NER", "NER", "NER", "NER", "NER", "NER", "NER", "BER", "MMR", "BER", "MMR", "MMR", "BER", "BER", "BER", "BER", "BER", "MMR", "BER", "POL", "MMR", "POL", "POL", "POL", "POL", "BER", "BER", "BER", "NER", "NER", "NER/BER")

genes <- unique(hmf_diplotypes$hgnc_symbol)

gene_pathway <- data.frame(hgnc_symbol = genes, pathway = collective_pathways, sub_pathway = specific_pathways)

`%notin%` <- Negate(`%in%`)


# Strict criteria

length(unique(hmf_diplotypes$sample))
deficient_strict_hmf_samples <- hmf_diplotypes[hmf_diplotypes$a1.max_score == 5 & hmf_diplotypes$a2.max_score == 5,]
length(unique(deficient_strict_hmf_samples$sample))
table(deficient_strict_hmf_samples$hgnc_symbol)

# adding the pathway information
deficient_strict_hmf_samples %<>% left_join(gene_pathway, by = "hgnc_symbol")
nrow(deficient_strict_hmf_samples[deficient_strict_hmf_samples$pathway == "MMR",])


# As you can see we will get only around 80 BER/NER deficient samples and it does not seem to be enough to do convincing statistics! For now I will just 
# follow the lenient criteria




# Lenient criteria

deficient_lenient_hmf_diplotypes <- hmf_diplotypes[hmf_diplotypes$biall_status == "deep_deletion" | hmf_diplotypes$diplotype_origin == "cnv_germ" & (hmf_diplotypes$a1.max_score==5 & hmf_diplotypes$a2.max_score>=4) | hmf_diplotypes$diplotype_origin == "cnv_som" & (hmf_diplotypes$a1.max_score==5 & hmf_diplotypes$a2.max_score>=3) | hmf_diplotypes$diplotype_origin %in% c('germ_som','som_som') & (hmf_diplotypes$a1.max_score==5 & hmf_diplotypes$a2.max_score==5), ]
# getting the samples with more than one driver mutation

dups <- deficient_lenient_hmf_diplotypes[duplicated(deficient_lenient_hmf_diplotypes$sample), "sample"]
single_mutant_deficient_lenient_hmf_diplotypes <- deficient_lenient_hmf_diplotypes[deficient_lenient_hmf_diplotypes$sample %notin% dups,]
single_mutant_deficient_lenient_hmf_diplotypes <- deficient_lenient_hmf_diplotypes
length(unique(deficient_lenient_hmf_diplotypes$sample))

table(deficient_lenient_hmf_diplotypes$hgnc_symbol)
# adding the pathway information
single_mutant_deficient_lenient_hmf_diplotypes %<>% left_join(gene_pathway, by = "hgnc_symbol")


# 12 samples were double mutants but both mutations were in the same pathway but 13 samples were double mutant meaning that had mutations in more than one pathway and therefore should be excluded!


# any(duplicated(single_mutant_deficient_lenient_hmf_diplotypes$sasample))
# in total 238 single mutant and 12 same pathway double mutant samples can be studied (250)
nrow(single_mutant_deficient_lenient_hmf_diplotypes)
nrow(double_same_mutant_deficient_lenient_hmf_diplotypes)/2
nrow(double_mutant_deficient_lenient_hmf_diplotypes)/2



# We also need to remove all the MMR samples as they are hypermutators! 67/238 and 9/12 samples are MMR deficient and need to be removed

sum(single_mutant_deficient_lenient_hmf_diplotypes$pathway == "MMR")
sum(double_same_mutant_deficient_lenient_hmf_diplotypes$pathway == "MMR")/2

# In general we obtained 171 mutant samples _ and 3 double-mutant samples 

# single_mutant_deficient_lenient_hmf_diplotypes <- single_mutant_deficient_lenient_hmf_diplotypes[single_mutant_deficient_lenient_hmf_diplotypes$pathway != "MMR",]


# removing the deficient samples to determine the grey listed samples!
hmf_diplotypes_without_deficients <- hmf_diplotypes[hmf_diplotypes$sample %notin% unique(single_mutant_deficient_lenient_hmf_diplotypes$sample),]
length(unique(hmf_diplotypes_without_deficients$sample))


## Grey list

geylist_lenient_hmf_samples <- hmf_diplotypes_without_deficients[hmf_diplotypes_without_deficients$diplotype_origin == "cnv_germ" & (hmf_diplotypes_without_deficients$a1.max_score==5 & hmf_diplotypes_without_deficients$a2.max_score==3) | hmf_diplotypes_without_deficients$diplotype_origin == "cnv_som" & (hmf_diplotypes_without_deficients$a1.max_score==5 & hmf_diplotypes_without_deficients$a2.max_score==2) | hmf_diplotypes_without_deficients$diplotype_origin %in% c('germ_som','som_som') & ((hmf_diplotypes_without_deficients$a1.max_score==5 & hmf_diplotypes_without_deficients$a2.max_score==4) | (hmf_diplotypes_without_deficients$a1.max_score==4 & hmf_diplotypes_without_deficients$a2.max_score==5) | (hmf_diplotypes_without_deficients$a1.max_score==4 & hmf_diplotypes_without_deficients$a2.max_score==4)), ]
length(unique(geylist_lenient_hmf_samples$sample))

# Proficient criteria



hmf_diplotypes_proficient <- hmf_diplotypes_without_deficients[hmf_diplotypes_without_deficients$sample %notin% unique(geylist_lenient_hmf_samples$sample),]
length(unique(hmf_diplotypes_proficient$sample))




##### Now I will annotate the metadata dataframe


metadata <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF/annotated_whitelisted_metadata.tsv",
                     sep = "\t", stringsAsFactors = F, header = T, )
metadata <- metadata[order(metadata$sampleId),]


metadata$pathway <- NA
metadata$sub_pathway <- NA
metadata$hgnc_symbol <- NA

metadata[metadata$sampleId %in% unique(hmf_diplotypes_proficient$sample), c("pathway", "sub_pathway")] <- "WT"
metadata[metadata$sampleId %in% unique(hmf_diplotypes_proficient$sample), c("hgnc_symbol")] <- "-"
metadata[metadata$sampleId %in% unique(geylist_lenient_hmf_samples$sample), c("pathway", "sub_pathway")] <- "Grey-listed"
metadata[metadata$sampleId %in% unique(geylist_lenient_hmf_samples$sample), c("hgnc_symbol")] <- "-"
metadata[metadata$sampleId %in% unique(single_mutant_deficient_lenient_hmf_diplotypes$sample), "pathway"] <- single_mutant_deficient_lenient_hmf_diplotypes[single_mutant_deficient_lenient_hmf_diplotypes$sample %in% metadata[metadata$sampleId %in% unique(single_mutant_deficient_lenient_hmf_diplotypes$sample), "sampleId"],"pathway"]
metadata[metadata$sampleId %in% unique(single_mutant_deficient_lenient_hmf_diplotypes$sample), "sub_pathway"] <- single_mutant_deficient_lenient_hmf_diplotypes[single_mutant_deficient_lenient_hmf_diplotypes$sample %in% metadata[metadata$sampleId %in% unique(single_mutant_deficient_lenient_hmf_diplotypes$sample), "sampleId"],"sub_pathway"]
metadata[metadata$sampleId %in% unique(single_mutant_deficient_lenient_hmf_diplotypes$sample), "hgnc_symbol"] <- single_mutant_deficient_lenient_hmf_diplotypes[single_mutant_deficient_lenient_hmf_diplotypes$sample %in% metadata[metadata$sampleId %in% unique(single_mutant_deficient_lenient_hmf_diplotypes$sample), "sampleId"],"hgnc_symbol"]

# These are samples that are MMR deicient or double mutants!
metadata[is.na(metadata$pathway), c("pathway", "sub_pathway")] <- "Grey-listed"
metadata[is.na(metadata$pathway), c("hgnc_symbol")] <- "-"


table(metadata$pathway)
table(metadata$sub_pathway)
table(metadata$hgnc_symbol)

# sum(is.na(metadata$pathway))


# write.table(metadata, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/labelled_metadata2.tsv", sep = "\t", quote = F, row.names = F)

metadata <- read.csv(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/labelled_metadata2.tsv", sep = "\t", header = T, stringsAsFactors = F)





# Overview of BER/NER/POL deficient samples across cancer types

## NER

png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/figs/NER-deficient.png")


xx <- table(metadata$primaryTumorLocation[metadata$pathway == "NER"])
tot_num <- sum(as.vector(xx))
names(xx) <- paste0(names(xx), "(", as.character(as.vector(xx)), "/", tot_num, ")")
pie(xx, cex = 0.75,  main = "Cancer-type distribution of NER-deficient samples (64 in total)")




dev.off()

# broken down to their gene level

png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/figs/NER-deficient_genes.png")


mm <- table(metadata$hgnc_symbol[metadata$pathway == "NER"])
tot_num <- sum(as.vector(mm))
names(mm) <- paste0(names(mm), "(", as.character(as.vector(mm)), "/42)")
colors <- RColorBrewer::brewer.pal(length(mm), "Spectral")
pie(mm, main = "Gene distribution of NER-deficient samples (42 in total)", col = colors)


dev.off()





# BER

png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/figs/BER-deficient.png")

xxx <- table(metadata$primaryTumorLocation[metadata$pathway == "BER"])
tot_num <- sum(as.vector(xxx))
names(xxx) <- paste0(names(xxx), "(", as.character(as.vector(xxx)), "/", tot_num, ")")
pie(xxx, cex = 0.75,  main = "Cancer-type distribution of BER-deficient samples (64 in total)")

dev.off()


# broken down to their gene level

png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/figs/BER-deficient_genes.png")


nn <- table(metadata$hgnc_symbol[metadata$pathway == "BER"])
tot_num <- sum(as.vector(nn))
names(nn) <- paste0(names(nn), "(", as.character(as.vector(nn)), "/", tot_num, ")")
yy <- sample(1:11, 11, replace = F)
colors <- RColorBrewer::brewer.pal(length(nn), "Spectral")
pie(nn[yy], clockwise = F, main = "Gene distribution of BER-deficient samples (74 in total)", col = colors)

dev.off()






## POL

png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/figs/POL-deficient.png")


xxxx <- table(metadata$primaryTumorLocation[metadata$pathway == "POL"])
tot_num <- sum(as.vector(xxxx))
names(xxxx) <- paste0(names(xxxx), "(", as.character(as.vector(xxxx)), "/", tot_num, ")")
pie(xxxx, cex = 0.75,  main = "Cancer-type distribution of POL-deficient samples (64 in total)")



dev.off()


# broken down to their gene level

png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/figs/POL-deficient_genes.png")


pp <- table(metadata$hgnc_symbol[metadata$pathway == "POL"])
tot_num <- sum(as.vector(pp))
names(pp) <- paste0(names(pp), "(", as.character(as.vector(pp)), "/", tot_num, ")")
colors <- RColorBrewer::brewer.pal(length(pp), "Spectral")
pie(pp, main = "Gene distribution of POL-deficient samples (48 in total)", col = colors)

dev.off()







## MMR

png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/figs/MMR-deficient.png")


xxxxx <- table(metadata$primaryTumorLocation[metadata$pathway == "MMR"])
tot_num <- sum(as.vector(xxxxx))
names(xxxxx) <- paste0(names(xxxxx), "(", as.character(as.vector(xxxxx)), "/", tot_num, ")")
pie(xxxxx, cex = 0.5, main = "Cancer-type distribution of MMR-deficient samples (64 in total)")


dev.off()


# broken down to their gene level

png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/ber-ner-status/figs/MMR-deficient_genes.png")


rr <- table(metadata$hgnc_symbol[metadata$pathway == "MMR"])
tot_num <- sum(as.vector(rr))
names(rr) <- paste0(names(rr), "(", as.character(as.vector(rr)), "/", tot_num, ")")
colors <- RColorBrewer::brewer.pal(length(rr), "Spectral")
pie(rr, main = "Gene distribution of MMR-deficient samples (64 in total)", col = colors)

dev.off()

