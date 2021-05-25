#!/bin/bash

VSEARCH=$(which vsearch)
THREADS=8
ENCODING=33

# list filenames containing R1

files=$(ls *-R1_001.fastq)

# loop over filenames and sort based on primer

for file_R1 in $files; do
	
	# get filename for R2
	file_R2=$(echo ${file_R1} | sed 's/R1/R2/')
	samp="${file_R1/-R1*/}"
	
	#write info to logfile
	echo "\n\nSample: ${samp} \n\n" >> ${samp}-merged.log

	# filter based on 18SV4 forward primer
  # Merge read pairs
  "${VSEARCH}" \
    --threads ${THREADS} \
    --fastq_mergepairs ${file_R1} \
    --reverse ${file_R2} \
    --fastq_ascii ${ENCODING} \
    --fastqout ${samp}-merged.fastq \
    --fastq_allowmergestagger \
    --quiet 2>> ${samp}-merged.log

done
