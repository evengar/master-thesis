#!/bin/bash

#SBATCH --account=nn9317k
## Job name:
#SBATCH --job-name=otu_contingency
## Number of tasks (aka processes) to start: Pure mpi, one cpu per task
#SBATCH --ntasks=1
## Amount of memory per cpu (= per task, since we get 1 cpu per task):
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1:00:00 
# turn on all mail notification, and also provide mail address:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evengar@student.ibv.uio.no
#SBATCH --output=%x-%u-%j.txt


set -o errexit #exit on errors

# remove all loaded modules
module purge
# OTU_contingency_table.py uses python 2.7
module load Python/2.7.18-GCCcore-10.2.0

FASTA="common-swarms/representatives_ns.fas"
SCRIPT="scripts/OTU_contingency_table.py"
STATS="${FASTA/.fas/.stats}"
SWARMS="${FASTA/.fas/.swarms}"
REPRESENTATIVES="${FASTA}"
UCHIME="common-swarms/common-cluster_1f_representatives.uchime"
ASSIGNMENTS="common-tax/common-ns-representatives.results"
QUALITY="common3/common.qual"
OTU_TABLE="common-tax/common-OTU.table"

python \
    "${SCRIPT}" \
    "${REPRESENTATIVES}" \
    "${STATS}" \
    "${SWARMS}" \
    "${UCHIME}" \
    "${QUALITY}" \
    "${ASSIGNMENTS}" \
    common-derep/*filtered.fasta > "${OTU_TABLE}"