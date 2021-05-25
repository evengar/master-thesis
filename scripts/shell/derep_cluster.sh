#!/bin/bash

#SBATCH --account=nn9317k
## Job name:
#SBATCH --job-name=derep_cluster
## Number of tasks (aka processes) to start: Pure mpi, one cpu per task
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
## Amount of memory per cpu (= per task, since we get 1 cpu per task):
#SBATCH --mem=8G
#SBATCH --time=08:00:00 
# turn on all mail notification, and also provide mail address:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evengar@student.ibv.uio.no
#SBATCH --output=%x-%u-%j.txt

set -o errexit #exit on errors

# remove all loaded modules
module purge
# load the vsearch module
module load VSEARCH/2.9.1-foss-2018b

VSEARCH=$(which vsearch)
SWARM=/cluster/projects/nn9317k/swarm/bin/swarm
TMP_FASTA=$(mktemp --tmpdir=".")
FINAL_FASTA="common-swarms/common-cluster.fas"

# Pool sequences
cat common-derep/*filtered.fasta > "${TMP_FASTA}"

# Dereplicate (vsearch)
"${VSEARCH}" --derep_fulllength "${TMP_FASTA}" \
             --sizein \
             --sizeout \
             --fasta_width 0 \
             --output "${FINAL_FASTA}" > /dev/null

rm -f "${TMP_FASTA}"

# Clustering
THREADS=1
TMP_REPRESENTATIVES=$(mktemp --tmpdir=".")
"${SWARM}" \
    -d 1 -f -t ${THREADS} -z \
    -i ${FINAL_FASTA/.fas/_1f.struct} \
    -s ${FINAL_FASTA/.fas/_1f.stats} \
    -w ${TMP_REPRESENTATIVES} \
    -o ${FINAL_FASTA/.fas/_1f.swarms} < ${FINAL_FASTA}

# Sort representatives
"${VSEARCH}" --fasta_width 0 \
             --sortbysize ${TMP_REPRESENTATIVES} \
             --output ${FINAL_FASTA/.fas/_1f_representatives.fas}
rm ${TMP_REPRESENTATIVES}
  
# Chimera checking
REPRESENTATIVES=${FINAL_FASTA/.fas/_1f_representatives.fas}
UCHIME=${REPRESENTATIVES/.fas/.uchime}
"${VSEARCH}" --uchime_denovo "${REPRESENTATIVES}" \
             --uchimeout "${UCHIME}"
