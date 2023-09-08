#!/bin/sh
#$ -S/bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ _M REPLACE_WITH_YOUR_INFO
#$ -pe whole_nodes 1
#########################################

module load python
module load bowtie
module load samtools
module load bedtools

#we source many variables from another shell script with needed variables.  Pass it when you call this script.
#do this by running sbatch process_seq.sh ./config.sh (or sub path to the config file if not in same folder).
source $1
raw_file=$2
f_name_start=${raw_file##*/}
bowtie_name=${Sorted_Save}${f_name_start/$Original_End/$Sorted_End}
processed_bowtie_name=${Sorted_Save}${f_name_start/$Original_End/$_processed.txt}

#first align, then compress to .bam files and sort aligned reads:
if [ ! -f "$bowtie_name" ]; then
	echo "Aligning ${f_name_start/$Original_End/} to $Genome"
	bowtie --threads 2 --sam $Bowtie_Args $Genome $raw_file | samtools view -u | samtools sort -T ${f_name_start/$Original_End/_sort_file.bam} -o $bowtie_name
fi

#if there is a mismatch at the 5' terminus, we trim that position (non-templated alignment)
#also discard all unmapped reads at this stage, and read of incorrect len (len < 14 or len > 45).
echo "Trimming the 5' end of all reads with 5' end mismatch"
srun process_bowtie_files.sh $bowtie_name $processed_bowtie_name

echo "making wig file:"
python ./make_wig.py $processed_bowtie_name
