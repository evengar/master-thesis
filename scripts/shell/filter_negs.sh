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

# pool amet negs
awk 'NR%2 == 0' *NC*amet* > neg_seqs_amet
# pool 18S negs
awk 'NR%2 == 0' *NC*18S* > neg_seqs_18S

for SAMP in S*amet*; do
  awk 'NR==FNR{id[$1]; next} !($0 in id)&&/^[^>]/{printf "%s\n%s\n", f, $0}{f=$0}' \
  neg_seqs_amet $SAMP > ${SAMP/fasta/filtered.fasta}
done

for SAMP in S*18S*; do
  awk 'NR==FNR{id[$1]; next} !($0 in id)&&/^[^>]/{printf "%s\n%s\n", f, $0}{f=$0}' \
  neg_seqs_18S $SAMP > ${SAMP/fasta/filtered.fasta}
done