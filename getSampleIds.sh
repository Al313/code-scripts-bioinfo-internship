#!/bin/bash



wd="/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj"
vcf_dir_hmf="/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics"
tadBedFileInput="$wd/annotation-beds/chromatin-status-tads.bed.bgz"  #args[3]
binBedFileInput="$wd/annotation-beds/chromatin-status-bins.bed.bgz"  #args[4]
repTimingFileInput="$wd/annotation-beds/repliseq_encode_mean_binned.bed.bgz"  #args[5]
repTimingFileInputBoxtel="$wd/annotation-beds/boxtel_paper_all_RepliSeq_median.bed.bgz" #args[6]
cpgInputFile="$wd/annotation-beds/cpgIslandExtUnmasked.bed.bgz"  #args[7]
repOrientationInputFile="$wd/annotation-beds/replication-direction-tableTerritories_Haradhvala_territories.rds.gz"  #args[8]

out_dir="$wd/HMF/annotated-variants"   #args[9]
job_dir="$out_dir/logs"
cohort="pole"

#samples.txt contains sample IDs obtained from spargling-genomics database. It currently contains only three sample IDsi!

sampleIds="/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/sample-ids/sampleIds.txt"
#counter=0

#numSamples=$(cat $wd/sample-ids/$sampleIds | wc -l)



touch sampleIdAvailable.txt
touch sampleIdNotAvailable.txt


for file_name in $(cat $sampleIds | grep -v ".[IV]$"); do

folder_name=$(echo ${file_name::-1})

path_to_vcf="${vcf_dir_hmf}/$(find ${vcf_dir_hmf} -maxdepth 1 | awk -v var=${folder_name} '$0~var' | rev | cut -d "/" -f 1 | rev | sort | head -n 1)/purple/${file_name}.purple.somatic.vcf.gz"

if [[ -f ${path_to_vcf} ]]; then
echo ${file_name} >> sampleIdAvailable.txt

else

echo ${file_name} >> sampleIdNotAvailable.txt

fi

done

