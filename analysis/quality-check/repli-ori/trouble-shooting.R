# +++ In this script I investigated how the replication orientation annotation is done for the variants +++

# Seemingly, the mut_strand from MutationalPatterns package annotates the variants for the replication strnd by taking the C/T conversion to A/G (like what it does for
# transcription strand annotation). Therefore we need to take that into account.


library(MutationalPatterns)
library(VariantAnnotation)
library(vcfR)



repli_strand_granges <- readRDS(file = "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/annotation-beds/replication-direction-tableTerritories_Haradhvala_territories.rds.gz")

path_to_vcf <- "/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics/170217_HMFregCPCT_FR12244678_FR12244592_CPCT02010350/purple/CPCT02010350T.purple.somatic.vcf.gz"


vcf <- read.vcfR(path_to_vcf)
df_SNVs <- cbind(as.data.frame(getFIX(vcf)), INFO2df(vcf))
df_SNVs <- df_SNVs[df_SNVs$FILTER == "PASS",]
df_SNVs<- df_SNVs[,c("CHROM", "POS", "REF", "ALT", "TNC", "PURPLE_AF", "KT", "MH", "REP_C", "REP_S")]
df_SNVs$POS <- as.integer(df_SNVs$POS)
chrOrder <-c((1:22),"X","Y")
df_SNVs <- df_SNVs[order(factor(df_SNVs$CHROM, levels = chrOrder, ordered = T), df_SNVs$POS),]

rownames(df_SNVs) <- paste0(df_SNVs$CHROM, ":", df_SNVs$POS, "_", df_SNVs$REF, "/", df_SNVs$ALT)
df_SNVs <- df_SNVs[isUnique(paste0(df_SNVs$CHROM,":", df_SNVs$POS)),]
colnames(df_SNVs)[c(1,2)] <- c("chromosome", "position")

# Remove indels!
df_SNVs <- df_SNVs[nchar(df_SNVs$ALT) == 1 & nchar(df_SNVs$REF) == 1,]

# A problem that came to my attention is that there are some variants call that point to exactly the same position in the genome but their PURPLE_AF position is different
# That can cause problem (specifically when running isOverlapping function). THerefore with the command below I get rid of them. For example in sample
# "CPCT02070491T" there were two rows for the same position that differ in ALT and PURPLE_AF fields (chrX:98875822)

#nrow(df_SNVs)
df_SNVs <- df_SNVs[isUnique(paste0(df_SNVs$chromosome,":",df_SNVs$position)),]

# First I let's sort the data frame based on the chromosome levels 
# (this step is needed for making GRanges and by doing sorting here the the order of rows will be the same for all the annotation)
chrOrder <-c((1:22),"X","Y")
df_SNVs <- df_SNVs[order(factor(df_SNVs$chromosome, levels = chrOrder, ordered = T), df_SNVs$position),]


# setting the row names
rownames(df_SNVs) <- paste0(df_SNVs$chromosome, ":", df_SNVs$position, "_", df_SNVs$REF, "/", df_SNVs$ALT)


# making a granges object
df_SNVs_gr <- GRanges(seqnames = df_SNVs$chromosome,ranges = IRanges(start = df_SNVs$position,end = df_SNVs$position),
                      ref = df_SNVs$REF, alt = df_SNVs$ALT) 
names(df_SNVs_gr) <- paste0(df_SNVs$chromosome, ":", df_SNVs$position, "_", df_SNVs$REF, "/", df_SNVs$ALT)



# setting the right seq levels
newStyle <- mapSeqlevels(seqlevels(df_SNVs_gr), "UCSC") 
df_SNVs_gr <- renameSeqlevels(df_SNVs_gr, newStyle)


rep_str_ann <- mut_strand(df_SNVs_gr, repli_strand_granges, mode = "replication")
df_SNVs_gr$repli_str_ori <- rep_str_ann


# Flip the replication strand orientation for variants that have G or A as their reference
index_for_right <- which((df_SNVs_gr$ref == "A" | df_SNVs_gr$ref == "G") & df_SNVs_gr$repli_str_ori == "right")

index_for_left <- which((df_SNVs_gr$ref == "A" | df_SNVs_gr$ref == "G") & df_SNVs_gr$repli_str_ori == "left")

df_SNVs_gr[index_for_right,]$repli_str_ori <- "left"
df_SNVs_gr[index_for_left,]$repli_str_ori <- "right"


# Regardless of what the REF nucleotide is, if a variant has an A or G as its ref nucleotide, the mut-_strand function flips its replication orientation strand:
df_SNVs_gr[df_SNVs_gr$ref == "A" & df_SNVs_gr$alt == "G"]
df_SNVs_gr[df_SNVs_gr$ref == "G" & df_SNVs_gr$alt == "A"]
df_SNVs_gr[df_SNVs_gr$ref == "T" & df_SNVs_gr$alt == "A"]
df_SNVs_gr[df_SNVs_gr$ref == "A" & df_SNVs_gr$alt == "T"]
df_SNVs_gr[df_SNVs_gr$ref == "T" & df_SNVs_gr$alt == "C"]
df_SNVs_gr[df_SNVs_gr$ref == "C" & df_SNVs_gr$alt == "T"][20:30]




# getting the annotated file
sample_id <- "CPCT02010350T"
vcf <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF/annotated-variants/pole/",
                              sample_id, "/annotated-", sample_id, ".txt"),
                header = T, sep = "\t", stringsAsFactors = F)

vcf[20:30,c("CHROM","POS","REF", "ALT", "rep_str_ann")]


#Inspecting the annotation in the reference bed file 
repli_strand_granges[seqnames(repli_strand_granges) == "chr1" & start(ranges(repli_strand_granges)) < 61552064 & end(ranges(repli_strand_granges)) > 61552064]






