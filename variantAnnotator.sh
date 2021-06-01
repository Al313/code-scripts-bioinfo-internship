#!/usr/bin/env bash


wd="/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/misc/processed/ali-proj"

tadBedFileInput="$wd/annotation-beds/chromatin-status-tads.bed.bgz"  #args[2]
binBedFileInput="$wd/annotation-beds/chromatin-status-bins.bed.bgz"  #args[3]
repTimingFileInput="$wd/annotation-beds/repliseq_encode_mean_binned.bed.bgz"  #args[4]
cpgInputFile="$wd/annotation-beds/cpgIslandExtUnmasked.bed.bgz"  #args[5]
repOrientationInputFile="$wd/annotation-beds/replication-direction-tableTerritories_Haradhvala_territories.rds.gz"  #args[6]

# out_dir="$wd/PCAWG/annotated-variants"
out_dir="$wd/HMF/annotated-variants" 

source ./.guix-profile-berner-proj/etc/profile

#counter=0
for i in $(cat $wd/sample-ids/samples.txt); do

	#if [ $counter -lt 1 ]; then

		if [[ -d $out_dir/$i ]]; then
			rm -rf $out_dir/$i;
			echo "Directory removed!";
		fi; 

		echo "Annotating $i";
		mkdir $out_dir/$i;
		Rscript $wd/variantAnnotator.R $i $tadBedFileInput $binBedFileInput $repTimingFileInput $cpgInputFile $repOrientationInputFile $out_dir;
	#	counter=$((${counter}+1));
	#fi;

done


