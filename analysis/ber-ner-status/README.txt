Link to diplotypes: /hpc/cuppen/projects/P0013_WGS_patterns_Diagn/datasets/processed/HMF_DR104_update3/scripts/gene_ann_ber_ner/diplotypes_germPonFilt.txt.gz

At the link below you can find all the infromation needed determining BER/NER status of samples:



This file contains the BER/NER status annotation for HMF and PCAWG samples:

/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/datasets/processed/HMF_DR104_update3/scripts/gene_ann_ber_ner/diplotypes_germPonFilt.txt.gz



Documentation can be found at:



Here you go: https://github.com/UMCUGenetics/cuplr/tree/master/geneDriverAnnotator

Luan:
"
Essentially, each gene is given a score out of 5 for allele 1 and allele 2 (a1, a2)
1=Benign, 5=Pathogenic
LOH (i.e. copy number loss of one allele) is always assigned to a1, and always has a score of 5 (edited) 
deep deletion (copy number loss of all alleles) is assigned to both a1 and a2, and a score of 5+5
"

Ali:
"
very nice! And just to be on the same page, what criteria do you consider for labelling a sample BER/NER deficient? e.g. both alleles => 4 then deficient, etc ... I can come up with one myself but I think it's better to use a reference criterion for later communication.
"

Luan:
"
The selection criteria are similar to what i used for training CHORD. Here's the script i used in the past: https://github.com/UMCUGenetics/CHORD/blob/master/training/sel_training_samples/sel_training_samples_HMF_DR010_DR047.R

biall_status=='deep_deletion' |
(diplotype_origin=='cnv_germ' & (a1.max_score==5 & a2.max_score>=4)) | ## Be more strict for germline SMNVs (on a2)
(diplotype_origin=='cnv_som' & (a1.max_score==5 & a2.max_score>=3)) | ## And less strict for somatic SMNVs (on a2)
(diplotype_origin %in% c('germ_som','som_som') & (a1.max_score==5 & a2.max_score==5)) ## Be very strict when dealing with 2x SMNV


Samples that don't meet the above criteria should be excluded (when analysing the target gene)
Also, samples that have MSI should be excluded (MSI is a hypermutator phenotype)

"




##### The sets of criteria I used for BER/NER labelling:



I used 4643 HMF samples! For cases that had more than one samples I used the first sample!

I marked samples that are POLE-mutant or associated with MSI and excluded them from the analysis. The reasoning for this exclusion is the fact that these samples are hypermutator and are associated with high TMB and might interfere with my analysis. 

No. of POLE-mutant samples: 31 (criterion: POLE be listed in driver genes file with a likelihood score of 0.5 or more) and 16 with lenient criterion I have shown below!
No. of MSI-high samples: 76


In addition I removed all samples that had more than 100,000 mutations (corresponding to 95th percentile of the snv counts of all data).
## Stats

mean =~ 9 Mut/Mb
median =~ 3.7 Mut/Mb


Set A (strict _ 1+1 groups): 


Deficient:

hmf_diplotypes$a1.max_score == 5 & hmf_diplotypes$a2.max_score == 5


Grey-listed:

hmf_diplotypes$a1.max_score => 4 $ hmf_diplotypes$a2.max_score == 4 | hmf_diplotypes$a1.max_score == 4 & hmf_diplotypes$a1.max_score => 4


Proficient:

Remaining samples!



## Stats

No. of BER-deficient:
No. of NER-deficient:
No. of grey-listed samples:




Set B (lenient _ 1+1 groups):


Deficient:

hmf_diplotypes$biall_status == "deep_deletion" | hmf_diplotypes$diplotype_origin == "cnv_germ" & (hmf_diplotypes$a1.max_score==5 & hmf_diplotypes$a2.max_score>=4) | hmf_diplotypes$diplotype_origin == "cnv_som" & (hmf_diplotypes$a1.max_score==5 & hmf_diplotypes$a2.max_score>=3) | hmf_diplotypes$diplotype_origin %in% c('germ_som','som_som') & (hmf_diplotypes$a1.max_score==5 & hmf_diplotypes$a2.max_score==5)

Grey-listed:

* Note that "hmf_diplotypes_without_deficients" is the metadata after the removal of deficient samples!

hmf_diplotypes_without_deficients$diplotype_origin == "cnv_germ" & (hmf_diplotypes_without_deficients$a1.max_score==5 & hmf_diplotypes_without_deficients$a2.max_score==3) | hmf_diplotypes_without_deficients$diplotype_origin == "cnv_som" & (hmf_diplotypes_without_deficients$a1.max_score==5 & hmf_diplotypes_without_deficients$a2.max_score==2) | hmf_diplotypes_without_deficients$diplotype_origin %in% c('germ_som','som_som') & ((hmf_diplotypes_without_deficients$a1.max_score==5 & hmf_diplotypes_without_deficients$a2.max_score==4) | (hmf_diplotypes_without_deficients$a1.max_score==4 & hmf_diplotypes_without_deficients$a2.max_score==5) | (hmf_diplotypes_without_deficients$a1.max_score==4 & hmf_diplotypes_without_deficients$a2.max_score==4)), ]


Proficient:

Remaining samples (MMR deficient samples were removed)!

## Stats

No. of BER-deficient: 74
No. of NER-deficient: 42
No. of grey-listed samples: 631
No. of POL-deficeint: 40


Set C (strict + 2+n groups):








Set D (lenient + 2+n groups):













** Interestingly enough, the POL genes are not detected as driver genes in the driver files of the hmf pipeline output. I checked the firs four POL labelled samples. Maybe not all POL genes are considered as drivers by the pipeline, but I know for sure that "POLE" is. However, when I checked the "CPCT02020762" sample which is labelled as POLE by our criteria, this gene was not listed in the driver file.
