#!/bin/sh
#########################################

#Change these variables to match what is needed in your project:

Genome="/bowtie_indices/?" #"," Seperated list of bowtie indexes to align to. Replace with your bowtie index
Raw_Data_Location="/raw_reads/" #Folder where the raw data lives.
Save_Folder_Location="/save_folder_loc/" #Folder for where to store all intermediary files.
Bowtie_Args="-m 1 -v 2"

########################################

Original_End=".fastq"
Sorted_Save=${Save_Folder_Location}aligned_and_sorted_files/
Sorted_End="_${Genome##*/}_sorted.bam"
Depth_Save=${Save_Folder_Location}depth_files/
Wig_Save=${Save_Folder_Location}wig_files/

#######################################

check_make_dir() {
	if [ ! -d "$1" ]; then
		echo "Making Directory $1"
		mkdir $1
	fi
}

check_make_dir "$Sorted_Save"
check_make_dir "$Depth_Save"
check_make_dir "$Wig_Save"
