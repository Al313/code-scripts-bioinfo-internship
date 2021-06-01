#!/bin/bash



wd="/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj"
tadBedFileInput="$wd/annotation-beds/chromatin-status-tads.bed.bgz"  #args[3]
binBedFileInput="$wd/annotation-beds/chromatin-status-bins.bed.bgz"  #args[4]
repTimingFileInput="$wd/annotation-beds/repliseq_encode_mean_binned.bed.bgz"  #args[5]
cpgInputFile="$wd/annotation-beds/cpgIslandExtUnmasked.bed.bgz"  #args[6]
repOrientationInputFile="$wd/annotation-beds/replication-direction-tableTerritories_Haradhvala_territories.rds.gz"  #args[7]

out_dir="$wd/HMF/annotated-variants"   #args[8]
job_dir="$out_dir/logs"


#samples.txt contains sample IDs obtained from spargling-genomics database. It currently contains only three sample IDs!

sampleIds="breast-ids.txt"
#counter=0

#numSamples=$(cat $wd/sample-ids/$sampleIds | wc -l)

for i in $(cat $wd/sample-ids/$sampleIds); do

path_to_vcf=

#if [ $counter -lt 3 ]; then

job_script=$job_dir/jobs/$i.sh


if [[ ! -p $job_script ]]; then
	touch $job_script;
fi


if [[ -d $out_dir/$i ]]; then
        rm -rf $out_dir/$i;
        echo "Directory removed!";
fi;

mkdir $out_dir/$i

echo "#!/bin/bash

#SBATCH --job-name=$(basename $job_script)
#SBATCH --output=$job_dir/outs/$i.txt
#SBATCH --error=$job_dir/errs/$i.txt
#SBATCH --time=00:30:00
#SBATCH --mem=10G

guixr load-profile $wd/.guix-profile-berner-proj --<<EOF

Rscript $wd/variantAnnotator.R $i $tadBedFileInput $binBedFileInput $repTimingFileInput $cpgInputFile $repOrientationInputFile $out_dir

EOF" > $job_script

if [[ ! -f ${job_script}.done ]]; then
	sbatch $job_script
else
	echo "Path for storing SBATCH outputs not available. Please construct the correct directory structure before running this program!"
fi

#counter=$((${counter}+1))

#else

#echo "$counter out of $numSamples samples selected"

#break

#fi

done








