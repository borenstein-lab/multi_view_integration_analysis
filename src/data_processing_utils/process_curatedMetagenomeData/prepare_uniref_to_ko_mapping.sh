#!/bin/bash

# The script reshapes humann's mapping from ko to UniRef90 from a wide format into a long format (with 2 columns: ko, UniRef90).
# Humann's mapping files are available here: https://github.com/biobakery/humann/tree/master/humann/data/misc

# Usage example: 
# cd /mnt/c/Users/efrat/Documents/GitHub/microbiome-project2/src/process_curatedMetagenomeData
# ./prepare_uniref_to_ko_mapping.sh /mnt/c/Users/efrat/Documents/REFERENCE_DBs/HUMANN2/utility_mapping/map_ko_uniref90.txt /mnt/c/Users/efrat/Documents/REFERENCE_DBs/HUMANN2/utility_mapping/map_ko_uniref90_long.txt

i_mapping_file_wide=$1
o_mapping_file_long=$2

while read -a LINE
do
	# Iterating over the array of UniRef codes, skipping the ko in the first position
	for ((i=1; i<${#LINE[@]}; i++));
	do
		printf  "${LINE[0]}\t${LINE[i]}\n"
	done
done < "$i_mapping_file_wide" > "$o_mapping_file_long"
