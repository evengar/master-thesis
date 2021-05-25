#!/bin/bash

#SBATCH --account=nn9317k
## Job name:
#SBATCH --job-name=cut_db
## Number of tasks (aka processes) to start: Pure mpi, one cpu per task
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=01:00:00 
# turn on all mail notification, and also provide mail address:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evengar@student.ibv.uio.no
#SBATCH --output=%x-%u-%j.txt

# Set options and load modules
set -o errexit #exit on errors

module purge
module load cutadapt/2.7-GCCcore-8.3.0-Python-3.7.4


# define objects
PRIMER_F=CGGTAAYTCCAGCTCYV # antimetazoa forward
PRIMER_R=TYRATCAAGAACGAAAGT # 18S reverse NB reverse complement

OUTPUT=PR2_13_ametF_18SR.fas
LOG="${OUTPUT/.fas/.log}"

MIN_LENGTH=32
MIN_F=$(( ${#PRIMER_F} * 2 / 3 ))	
MIN_R=$(( ${#PRIMER_R} * 2 / 3 ))

CUTADAPT="$(which cutadapt) --discard-untrimmed --minimum-length ${MIN_LENGTH}"

dos2unix < pr2_version_4.13.0_18S_UTAX.fasta | \
    sed '/^>/ s/;tax=k:/ /
         /^>/ s/,[dpcofgs]:/|/g
         /^>/ ! s/U/T/g' | \
    ${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2> "${LOG}" | \
    ${CUTADAPT} -a "${PRIMER_R}" -O "${MIN_R}" - 2>> "${LOG}" > "${OUTPUT}"