#!/bin/bash

#SBATCH --account=nn9317k
## Job name:
#SBATCH --job-name=local_derep
## Number of tasks (aka processes) to start: Pure mpi, one cpu per task
#SBATCH --ntasks=1
## Amount of memory per cpu (= per task, since we get 1 cpu per task):
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00 
# turn on all mail notification, and also provide mail address:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evengar@student.ibv.uio.no
#SBATCH --output=%x-%u-%j.txt

INPUT=common-cluster_1f_representatives.fas
OUTPUT=representatives_ns.fas

awk 'BEGIN{FS="[=;]"}NR%2==1&&$3>1{print;getline;print}' $INPUT > $OUTPUT

awk 'BEGIN{FS="[;]"}NR%2==1 {print substr($1, 2)}' $OUTPUT > seeds_filtered

awk 'NR==FNR{id[$1]; next} $3 in id' seeds_filtered ${INPUT/_representatives.fas/.stats} > ${OUTPUT/fas/stats}

awk 'BEGIN{FS=";"} NR==FNR{id[$1]; next} $1 in id' seeds_filtered ${INPUT/_representatives.fas/.swarms} > ${OUTPUT/fas/swarms}

awk 'NR==FNR{id[$1]; next} $1 in id' seeds_filtered ${INPUT/_representatives.fas/.struct} > ${OUTPUT/fas/struct}


