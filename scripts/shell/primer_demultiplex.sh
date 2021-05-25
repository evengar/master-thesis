#!/bin/bash

#SBATCH --account=nn9317k
## Job name:
#SBATCH --job-name=primer_demultiplex
## Number of tasks (aka processes) to start: Pure mpi, one cpu per task
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
## Amount of memory per cpu (= per task, since we get 1 cpu per task):
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00 
# turn on all mail notification, and also provide mail address:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evengar@student.ibv.uio.no
#SBATCH --output=%x-%u-%j.txt

# Set options and load modules
set -o errexit #exit on errors

module purge
module load cutadapt/2.7-GCCcore-8.3.0-Python-3.7.4

JOB=common3/primer_demultiplex_new
echo "" >> ${JOB}.log

	
p18S_R=TYRATCAAGAACGAAAGT
amet_F=CGGTAAYTCCAGCTCYV
amet_R=CGCAAGRSTGAAAYTTAAAG

# specify minimum overlap
MIN_amet_F=$(( ${#amet_F} * 2 / 3 ))	
MIN_amet_R=$(( ${#amet_R} * 2 / 3 ))
MIN_p18S_R=$(( ${#p18S_R} * 2 / 3 ))

# make temp objects for storing intermediary files
TMP_18S=$(mktemp)
TMP_amet=$(mktemp)
TMP_all=$(mktemp)

# specify number of cores (this can be higher on Saga)
CORES=4

# list files
files=$(ls merged-fastq/*-merged.fastq)

# loop over filenames and sort based on primer

for file in $files; do
	# get sample name
	samp="$(echo $file | sed -r 's/^.*\/(.*)-merged.fastq/\1/')"
	# write to logfile
	echo "Sample ${samp}:" >> ${JOB}.log
	
	# cutting by amet_F first
	cutadapt -g $amet_F \
	--discard-untrimmed \
	--cores $CORES \
	-O $MIN_amet_F \
	$file > $TMP_all 2>> $JOB.log
	
	
	# demultiplexing by amet_R-primer
	cutadapt -a $amet_R \
	--untrimmed-output $TMP_18S \
	--cores $CORES \
	-O $MIN_amet_R \
	$TMP_all > $TMP_amet 2>> $JOB.log
	
	# cut 18S
	cutadapt -a $p18S_R \
	--discard-untrimmed \
	--minimum-length 300 \
	--maximum-length 450 \
	--cores $CORES \
	-O $MIN_p18S_R \
	$TMP_18S > common3/${samp}-18S-merged.fastq 2>> ${JOB}.log
	
	# cut amet
	cutadapt -a $p18S_R \
	--cores $CORES \
	--minimum-length 300 \
	--maximum-length 450 \
	--discard-untrimmed \
	-O $MIN_p18S_R \
	$TMP_amet > common3/${samp}-amet-merged.fastq 2>> ${JOB}.log

done

# discard temp files
rm -f $TMP_18S $TMP_amet
