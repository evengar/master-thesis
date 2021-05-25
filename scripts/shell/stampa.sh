#!/bin/bash

#SBATCH --account=nn9317k
## Job name:
#SBATCH --job-name=tax_assign
## Number of tasks (aka processes) to start: Pure mpi, one cpu per task
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
## Amount of memory per cpu (= per task, since we get 1 cpu per task):
#SBATCH --mem=4G
#SBATCH --time=3:00:00 
# turn on all mail notification, and also provide mail address:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evengar@student.ibv.uio.no
#SBATCH --output=%x-%u-%j.txt


set -o errexit #exit on errors

# remove all loaded modules
module purge
# load the vsearch module
module load VSEARCH/2.9.1-foss-2018b

# variables
QUERY="common-swarms/representatives_ns.fas"
DATABASE="pr2/PR2_13_ametF_18SR.fas"
VSEARCH=$(which vsearch)


# search for best hits
$VSEARCH \
    --usearch_global ${QUERY} \
    --dbmask none \
    --qmask none \
    --rowlen 0 \
    --notrunclabels \
    --userfields query+id1+target \
    --maxaccepts 0 \
    --maxrejects 32 \
    --top_hits_only \
    --output_no_hits \
    --db ${DATABASE} \
    --id 0.5 \
    --iddef 1 \
    --userout - | sed 's/;size=/_/ ; s/;//' > common-tax/hits.representatives

module purge
# stampa_merge.py uses python 3
module load Python/3.8.6-GCCcore-10.2.0

# in case of multi-best hit, find the last-common ancestor
python scripts/stampa_merge.py $(pwd)/common-tax

# sort by decreasing abundance
sort -k2,2nr -k1,1d common-tax/results.representatives > common-tax/common-ns-representatives.results
