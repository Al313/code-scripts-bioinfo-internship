############ Automated script

library(dplyr)
library(magrittr)
library(ggplot2)
library(stringr)

sample_ids <- scan("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/annotated-pole-sample-ids.txt",
                   what = character()) # These are sample ids of pan-cancer pole mutants

sample_ids <- scan("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/colorectum-pole-sample-ids.txt",
                   what = character()) # These are sample ids of colorectal pole mutants

sample_ids <- scan("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/sample-ids/cancer-type/Lung/annotated.txt",
                   what = character()) # These are sample ids of lung cancers that I want to use as negative control for pole asymmetry
sample_ids <- sample_ids[!str_detect(sample_ids, pattern = "^21")]


myDataFrame <- data.frame("c>a_unknown" = NA, "c>a_left" = NA, "c>a_right" = NA, "g>t_unknown" = NA, "g>t_left" = NA, "g>t_right" = NA,
                          "c>g_unknown" = NA, "c>g_left" = NA, "c>g_right" = NA, "g>c_unknown" = NA, "g>c_left" = NA, "g>c_right" = NA,
                          "c>t_unknown" = NA, "c>t_left" = NA, "c>t_right" = NA, "g>a_unknown" = NA, "g>a_left" = NA, "g>a_right" = NA,
                          "t>a_unknown" = NA, "t>a_left" = NA, "t>a_right" = NA, "a>t_unknown" = NA, "a>t_left" = NA, "a>t_right" = NA,
                          "t>c_unknown" = NA, "t>c_left" = NA, "t>c_right" = NA, "a>g_unknown" = NA, "a>g_left" = NA, "a>g_right" = NA,
                          "t>g_unknown" = NA, "t>g_left" = NA, "t>g_right" = NA, "a>c_unknown" = NA, "a>c_left" = NA, "a>c_right" = NA,
                          "tot_unknown" = NA, "tot_left" = NA, "tot_right" = NA)


for (i in 1:length(sample_ids)){
  print(i)
  sample_id <- sample_ids[i]
  vcf <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF/annotated-variants/cancer-type/Lung/", sample_id,
                                "/annotated-", sample_id, ".txt"),
                  header = T, sep = "\t", stringsAsFactors = F)
  
  vcf1 <- vcf[vcf$REF == "C" & vcf$ALT == "A", ]
  table1 <- table(vcf1$rep_str_ann)
  myDataFrame[i,1:3] <- as.numeric(table1)
  
  vcf2 <- vcf[vcf$REF == "G" & vcf$ALT == "T", ]
  table2 <- table(vcf2$rep_str_ann)
  myDataFrame[i,4:6] <- as.numeric(table2)
  
  vcf3 <- vcf[vcf$REF == "C" & vcf$ALT == "G", ]
  table3 <- table(vcf3$rep_str_ann)
  myDataFrame[i,7:9] <- as.numeric(table3)
  
  vcf4 <- vcf[vcf$REF == "G" & vcf$ALT == "C", ]
  table4 <- table(vcf4$rep_str_ann)
  myDataFrame[i,10:12] <- as.numeric(table4)
  
  vcf5 <- vcf[vcf$REF == "C" & vcf$ALT == "T", ]
  table5 <- table(vcf5$rep_str_ann)
  myDataFrame[i,13:15] <- as.numeric(table5)
  
  vcf6 <- vcf[vcf$REF == "G" & vcf$ALT == "A", ]
  table6 <- table(vcf6$rep_str_ann)
  myDataFrame[i,16:18] <- as.numeric(table6)
  
  vcf7 <- vcf[vcf$REF == "T" & vcf$ALT == "A", ]
  table7 <- table(vcf7$rep_str_ann)
  myDataFrame[i,19:21] <- as.numeric(table7)
  
  vcf8 <- vcf[vcf$REF == "A" & vcf$ALT == "T", ]
  table8 <- table(vcf8$rep_str_ann)
  myDataFrame[i,22:24] <- as.numeric(table8)
  
  vcf9 <- vcf[vcf$REF == "T" & vcf$ALT == "C", ]
  table9 <- table(vcf9$rep_str_ann)
  myDataFrame[i,25:27] <- as.numeric(table9)
  
  vcf10 <- vcf[vcf$REF == "A" & vcf$ALT == "G", ]
  table10 <- table(vcf10$rep_str_ann)
  myDataFrame[i,28:30] <- as.numeric(table10)
  
  vcf11 <- vcf[vcf$REF == "T" & vcf$ALT == "G", ]
  table11 <- table(vcf11$rep_str_ann)
  myDataFrame[i,31:33] <- as.numeric(table11)
  
  vcf12 <- vcf[vcf$REF == "A" & vcf$ALT == "C", ]
  table12 <- table(vcf12$rep_str_ann)
  myDataFrame[i,34:36] <- as.numeric(table12)
  
  table13 <- table(vcf$rep_str_ann)
  myDataFrame[i,37:39] <- as.numeric(table13)
  
  rownames(myDataFrame[i,]) <- sample_id
  
}

# getting the total number of mutations that are annotated for their replication orientation that can be used for correction

myDataFrame %<>% mutate(tot_c.a = c.a_left + c.a_right)
myDataFrame %<>% mutate(tot_g.t = g.t_left + g.t_right)
myDataFrame %<>% mutate(tot_c.g = c.g_left + c.g_right)
myDataFrame %<>% mutate(tot_g.c = g.c_left + g.c_right)
myDataFrame %<>% mutate(tot_c.t = c.t_left + c.t_right)
myDataFrame %<>% mutate(tot_g.a = g.a_left + g.a_right)
myDataFrame %<>% mutate(tot_t.a = t.a_left + t.a_right)
myDataFrame %<>% mutate(tot_a.t = a.t_left + a.t_right)
myDataFrame %<>% mutate(tot_t.c = t.c_left + t.c_right)
myDataFrame %<>% mutate(tot_a.g = a.g_left + a.g_right)
myDataFrame %<>% mutate(tot_t.g = t.g_left + t.g_right)
myDataFrame %<>% mutate(tot_a.c = a.c_left + a.c_right)
myDataFrame %<>% mutate(tot_all = tot_left + tot_right)


# correcting for number of mutations
myDataFrame[,c(2,3)] <- myDataFrame[,c(2,3)]/myDataFrame[,40]
myDataFrame[,c(5,6)] <- myDataFrame[,c(5,6)]/myDataFrame[,41]
myDataFrame[,c(8,9)] <- myDataFrame[,c(8,9)]/myDataFrame[,42]
myDataFrame[,c(11,12)] <- myDataFrame[,c(11,12)]/myDataFrame[,43]
myDataFrame[,c(14,15)] <- myDataFrame[,c(14,15)]/myDataFrame[,44]
myDataFrame[,c(17,18)] <- myDataFrame[,c(17,18)]/myDataFrame[,45]
myDataFrame[,c(20,21)] <- myDataFrame[,c(20,21)]/myDataFrame[,46]
myDataFrame[,c(23,24)] <- myDataFrame[,c(23,24)]/myDataFrame[,47]
myDataFrame[,c(26,27)] <- myDataFrame[,c(26,27)]/myDataFrame[,48]
myDataFrame[,c(29,30)] <- myDataFrame[,c(29,30)]/myDataFrame[,49]
myDataFrame[,c(32,33)] <- myDataFrame[,c(32,33)]/myDataFrame[,50]
myDataFrame[,c(35,36)] <- myDataFrame[,c(35,36)]/myDataFrame[,51]

myDataFrame[,c(38,39)] <- myDataFrame[,c(38,39)]/myDataFrame[,52]



plotdf <- data.frame("mutation_type" = NA, "replication_orientation" = NA, "log2ratio" = NA)
counter <- 0

for (i in c("a","c")) {
  
  i_holder <- i
  
  for (j in c("a","c","t","g")) {
    
    i <- i_holder
    
    if (i != j) {
      j_holder <- j
      for (k in c("left", "right")) {
        i <- i_holder
        j <- j_holder
        txt1 <- paste0(i, ".", j)
        # print(mean(myDataFrame[,paste0(i, ".", j, "_left")]))
        
        # i complementary
        
        if (i == "a"){
          i = "t"
        } else if (i == "t") {
          i = "a"
        } else if (i == "c") {
          i = "g"
        } else if (i == "g") {
          i = "c"
        }
        
        
        # j complementary
        if (j == "a"){
          j = "t"
        } else if (j == "t") {
          j = "a"
        } else if (j == "c") {
          j = "g"
        } else if (j == "g") {
          j = "c"
        }
        
        txt2 <- paste0(i, ".", j)
        
        log2ratio <- log2(mean(myDataFrame[,paste0(txt1, "_", k)])/mean(myDataFrame[,paste0(txt2, "_", k)]))
        
        counter <- counter + 1
        plotdf[counter,] <- c(paste0(txt1,"/", txt2), k, log2ratio)
        
        print(paste0("log2 ratio of ", txt1, "_", k, " and ", txt2, "_", k, " is ", log2ratio))
        # print(mean(myDataFrame[,paste0(i, ".", j, "_left")]))
        
      }
      
    }
  }
}


########

# Pan-cancer pole mutants

# saveRDS(myDataFrame, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/pole-mutation-count.rds")
myDataFrame <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/pole-mutation-count.rds")

myDataFrame_cp <- myDataFrame
myDataFrame <- myDataFrame_cp


# Visualiation


# png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/boxplot-rep-strand-bias-all-pole.png")
# 
# boxplot(myDataFrame[,c(2,3,5,6)], main = "POLE samples's replication strand bias", ylab = "observed/expected replication bias", xlab = "Mutation type")
# dev.off()
# 
# 
# png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/barplot-rep-strand-bias-all-pole.png")
# 
# sdd <- apply(myDataFrame[,c(2,3,5,6)], 2, sd, na.rm = T)/sqrt(nrow(myDataFrame[,c(2,3,5,6)]))
# summ <- apply(myDataFrame[,c(2,3,5,6)], 2, mean, na.rm = T)
# 
# mid <- barplot(summ, col = c("red", "orange", "lightblue", "darkblue"), width = 0.4, space = 0.75, main = "Replication strand bias in POLE Samples",
#                ylab = "observed/expected replication bias", xlab = "Mutation types")
# 
# arrows(x0=mid, y0=summ-sdd, x1=mid, y1=summ+sdd, code=3, col="black", lwd=2, angle=90, length=0.05)
# dev.off()


# Right-replicating genome
png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/right-replicating-barplot-rep-strand-bias-all-pole.png")
values <- as.numeric(plotdf[plotdf$replication_orientation == "right",c(3)])
names(values) <- toupper(gsub("[.]", ">", plotdf[plotdf$replication_orientation == "right",c(1)]))


barplot(values, ylim = c(-1,1), col = rainbow(6), main = "Right-replicating regions for pan-cancer POLE-mutants",
        ylab = "log2 ratio")
dev.off()


# Left-replicating genome
png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/left-replicating-barplot-rep-strand-bias-all-pole.png")

values <- as.numeric(plotdf[plotdf$replication_orientation == "left",c(3)])
names(values) <- toupper(gsub("[.]", ">", plotdf[plotdf$replication_orientation == "left",1]))


barplot(values, ylim = c(-1,1), col = rainbow(6), main = "Left-replicating regions for pan-cancer POLE-mutants",
        ylab = "log2 ratio")

dev.off()


# overall
png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/entire-genome-barplot-rep-strand-bias-all-pole.png")
values <- as.numeric(plotdf[plotdf$replication_orientation == "right",c(3)]) + as.numeric(plotdf[plotdf$replication_orientation == "left",c(3)])
names(values) <- toupper(gsub("[.]", ">", plotdf[plotdf$replication_orientation == "left",c(1)]))

barplot(values, ylim = c(-1,1), col = rainbow(6), main = "Entire genome for pan-cancer POLE-mutants",
        ylab = "log2 ratio")
dev.off()


########

# pole-colorectum-sample-ids.txt

# saveRDS(myDataFrame, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/pole-colorectum-mutation-count.rds")
myDataFrame <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/pole-colorectum-mutation-count.rds")

myDataFrame_cp2 <- myDataFrame
myDataFrame <- myDataFrame_cp2



# visualization

# png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/boxplot-rep-strand-bias-colorectum.png")
# boxplot(myDataFrame[,c(2,3,5,6)], main = "POLE samples's replication strand bias (Colorectal Cohort)", ylab = "observed/expected replication bias", xlab = "Mutation type")
# dev.off()
# 
# sdd <- apply(myDataFrame[,c(2,3,5,6)], 2, sd, na.rm = T)/sqrt(nrow(myDataFrame[,c(2,3,5,6)]))
# summ <- apply(myDataFrame[,c(2,3,5,6)], 2, mean, na.rm = T)
# 
# png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/barplot-trp-strand-bias-colorectum.png")
# 
# mid <- barplot(summ, col = c("red", "orange", "lightblue", "darkblue"), width = 0.4, space = 0.75, main = "Replication Strand Bias in POLE Samples (Colorectal Cohort)",
#                ylab = "Observed/expected replication bias", xlab = "Mutation types")
# arrows(x0=mid, y0=summ-sdd, x1=mid, y1=summ+sdd, code=3, col="black", lwd=2, angle=90, length=0.05)
# dev.off()



png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/right-replicating-barplot-rep-strand-bias-colorectum.png")
values <- as.numeric(plotdf[plotdf$replication_orientation == "right",c(3)])
names(values) <- toupper(gsub("[.]", ">", plotdf[plotdf$replication_orientation == "right",c(1)]))


barplot(values, ylim = c(-1,1), col = rainbow(6), main = "Right-replicating regions for colorectum POLE-mutants",
        ylab = "log2 ratio")
dev.off()



png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/left-replicating-barplot-rep-strand-bias-colorectum.png")

values <- as.numeric(plotdf[plotdf$replication_orientation == "left",c(3)])
names(values) <- toupper(gsub("[.]", ">", plotdf[plotdf$replication_orientation == "left",c(1)]))

barplot(values, ylim = c(-1,1), col = rainbow(6), main = "Left-replicating regions for colorectum POLE-mutants",
        ylab = "log2 ratio")

dev.off()


# overall
png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/entire-genome-barplot-rep-strand-bias-colorectum.png")
values <- as.numeric(plotdf[plotdf$replication_orientation == "right",c(3)]) + as.numeric(plotdf[plotdf$replication_orientation == "left",c(3)])
names(values) <- toupper(gsub("[.]", ">", plotdf[plotdf$replication_orientation == "left",c(1)]))

barplot(values, ylim = c(-1,1), col = rainbow(6), main = "Entire genome for colorectum POLE-mutants",
        ylab = "log2 ratio")
dev.off()







########

# lung-cancer not pole mutant


# saveRDS(myDataFrame, file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/non-pole-lung.rds")
myDataFrame <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/non-pole-lung.rds")



png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/right-replicating-barplot-rep-strand-bias-non-pole-lung.png")
values <- as.numeric(plotdf[plotdf$replication_orientation == "right",c(3)])
names(values) <- toupper(gsub("[.]", ">", plotdf[plotdf$replication_orientation == "right",c(1)]))


barplot(values, ylim = c(-1,1), col = rainbow(6), main = "Right replicating regions for non-POLE lung cancers",
        ylab = "log2 ratio")
dev.off()

png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/left-replicating-barplot-rep-strand-bias-non-pole-lung.png")

values <- as.numeric(plotdf[plotdf$replication_orientation == "left",c(3)])
names(values) <- toupper(gsub("[.]", ">", plotdf[plotdf$replication_orientation == "left",c(1)]))

barplot(values, ylim = c(-1,1), col = rainbow(6), main = "Left replicating regions for non-POLE lung cancers",
        ylab = "log2 ratio")

dev.off()


png(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori/entire-genome-barplot-rep-strand-bias-non-pole-lung.png")
values <- as.numeric(plotdf[plotdf$replication_orientation == "right",c(3)]) + as.numeric(plotdf[plotdf$replication_orientation == "left",c(3)])
names(values) <- toupper(gsub("[.]", ">", plotdf[plotdf$replication_orientation == "left",c(1)]))

barplot(values, ylim = c(-1,1), col = rainbow(6), main = "Entire genome for pan-cancer POLE-mutants",
        ylab = "log2 ratio")
dev.off()

