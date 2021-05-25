#!/bin/bash

#SBATCH --account=nn9317k
## Job name:
#SBATCH --job-name=local_derep
## Number of tasks (aka processes) to start: Pure mpi, one cpu per task
#SBATCH --ntasks=1
## Amount of memory per cpu (= per task, since we get 1 cpu per task):
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00 
# turn on all mail notification, and also provide mail address:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evengar@student.ibv.uio.no
#SBATCH --output=%x-%u-%j.txt

set -o errexit #exit on errors

# remove all loaded modules
module purge
# load the vsearch module
module load VSEARCH/2.9.1-foss-2018b

# Set variables
FILES=$(ls common3/*-merged.fastq)
TMP_FASTQ=$(mktemp)
TMP_FASTA=$(mktemp)
OUTPUT=$(mktemp)
LOG="local_derep_common.log"
VSEARCH=$(which vsearch)
QUALITY_FILE="common3/common.qual"


for FILE in $FILES; do

	SAMPLE=$(echo $FILE | sed -r 's/^.*\/(.*)-merged.fastq/\1/')
	FINAL_FASTA="common-derep/${SAMPLE}-derep.fasta"


	# Discard sequences containing Ns, add expected error rates
	# This is for assembling the quality file later
	"$VSEARCH" \
		--quiet \
		--fastq_filter "${FILE}" \
		--fastq_maxns 0 \
		--relabel_sha1 \
		--eeout \
		--fastqout "${TMP_FASTQ}" 2>> "${LOG}"

	# Discard sequences containing Ns, convert to fasta
	"${VSEARCH}" \
		--quiet \
		--fastq_filter "${FILE}" \
		--fastq_maxns 0 \
		--fastaout "${TMP_FASTA}" 2>> "${LOG}"

	# Dereplicate at the study level
	"${VSEARCH}" \
		--quiet \
		--derep_fulllength "${TMP_FASTA}" \
		--sizeout \
		--fasta_width 0 \
		--relabel_sha1 \
		--output "${FINAL_FASTA}" 2>> "${LOG}"
		
	# Discard quality lines, extract hash, expected error rates and read length
	sed 'n;n;N;d' "${TMP_FASTQ}" | \
		awk 'BEGIN {FS = "[;=]"}
			 {if (/^@/) {printf "%s\t%s\t", $1, $3} else {print length($1)}}' | \
		tr -d "@" >> "${OUTPUT}"
done
 
# Produce the final quality file
sort -k3,3n -k1,1d -k2,2n "${OUTPUT}" | \
    uniq --check-chars=40 > "${QUALITY_FILE}"

# Clean
rm -f "${TMP_FASTQ}" "${TMP_FASTA}" "${OUTPUT}"