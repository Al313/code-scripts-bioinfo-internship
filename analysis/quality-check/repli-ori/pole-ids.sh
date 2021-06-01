#!/bin/bash


metadata="/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF/metadata.tsv"
wd="/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/analysis/quality-check/repli-ori"
hmf="/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics"

#for sample_id in $(cat $wd/sample-ids.txt | grep -v ".[IV]$"); do echo ${sample_id} >> ${wd}/sample-ids-one-per-patient.txt; done



for sample_id in $(cat ${wd}/colorectum-sample-ids-one-per-patient.txt); do

folder_name=$(echo ${sample_id::-1})

path_to_drivers="${hmf}/$(find ${hmf} -maxdepth 1 | awk -v var=${folder_name} '$0~var' | rev | cut -d "/" -f 1 | rev | sort | head -n 1)/linx/${sample_id}.driver.catalog.tsv"

if [[ $(cat ${path_to_drivers} | grep "POLE" | wc -l) -gt 0 ]]; then if (( $(echo "$(cat ${path_to_drivers} | grep "POLE" | cut -f 7) > 0.5" | bc -l) )); then echo ${sample_id} >> ${wd}/colorectum-pole-sample-ids.txt; fi; fi


done


