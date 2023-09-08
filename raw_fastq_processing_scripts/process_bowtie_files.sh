#!/bin/sh
#$ -S/bin/sh
#$ -cwd
#$ -V
#$ -m e
#$ _M YOUR_INFO_HERE
#$ -pe whole_nodes 1
#########################################

module load samtools
#Pass the name of the bowtie file we are aligning to this when running.
#Run this script using: sbatch process_bowtie_files.sh path_to_bowtie_file
file=$1
new_file=$2
echo "Processing Aligned file: $file \n Discarding reads with len < 14 or len >45. \n If first read doesn't align due to potential non-templated addition, trim that read. \n Prooduce compressed info file. "
samtools view $file -F 0x4 | awk 'BEGIN{
	bases = "A T G C"
	split(bases, temp)
	for (i in temp)
		bases_list[temp[i]]
} {
	len_seq=length($10)
	if ((len_seq > 14 && len_seq < 45) && !((substr($(NF-2), 6, 1) == '0' && (substr($(NF-2), 8, 1) in bases_list))))
		five_end=$4
		strand="+"
		chrom=$3
		if ($2 == "16")
			strand="-"
		if (substr($(NF-2), 6, 1) == '0'){
			if (strand == "+")
				five_end = five_end + 1
			len_seq = len_seq-1
		}
		other_end = five_end + len_seq
		if (strand == "+")
			print $3 "\t" strand "\t" five_end "\t" other_end
		else
			print $3 "\t" strand "\t" other_end "\t" five_end

}
' > $new_file
echo "done!"
