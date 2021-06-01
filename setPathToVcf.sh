
wd="/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj/HMF"
inPutDir="/hpc/cuppen/shared_resources/HMF_data/DR-104-update3/somatics"

cd $wd


rowNo=$(cat metadata_whitelisted.tsv | sed -e '1d' | wc -l)

primaryTumorLocs=$(cat metadata_whitelisted.tsv | sed -e '1d' | cut -d$'\t' -f 7)
sampleNames=$(cat metadata_whitelisted.tsv | sed -e '1d' | cut -d$'\t' -f 3)
sampleIds=$(cat metadata_whitelisted.tsv | sed -e '1d' | cut -d$'\t' -f 2)



for i in $(eval echo "{1..$rowNo}"); do

primaryTumorLoc=$(echo ${primaryTumorLocs} | cut -d " " -f ${i})
sampleName=$(echo ${sampleNames} | cut -d " " -f ${i})
sampleId=$(echo ${sampleIds} | cut -d " " -f ${i})



path_to_vcf="${inPutDir}/${sampleName}/purple/${sampleId}.purple.somatic.vcf.gz"

if [ ! -f ${path_to_vcf} ]; then
    echo "${path_to_vcf} does not exist."
fi

echo $i
done
