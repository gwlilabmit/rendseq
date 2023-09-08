#!/bin/sh
#$ -S/bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ _M REPLACE_WITH_YOUR_INFO
#$ -pe whole_nodes 1
#########################################

#We source many variables from another shell script with needed variables.  Pass it when you call this script.
#Run this script using: sbatch align_all_files.sh path_to_config_file
source $1

for f in ${Raw_Data_Location}*$Original_End; do
	echo "Starting Script for $f ..."
	sbatch raw_fastq_align.sh $1 $f
done
