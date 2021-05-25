#!/bin/bash

#SBATCH --account=nn9317k
#SBATCH --job-name=self_match
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=4G
#SBATCH --time=2:00:00 
# turn on all mail notification, and also provide mail address:
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evengar@student.ibv.uio.no
#SBATCH --output=%x-%u-%j.txt

module purge
module load VSEARCH/2.9.1-foss-2018b

REPRESENTATIVES=common-swarms/representatives_ns.fas
OUTPUT=lulu/match_list_common.txt

vsearch --usearch_global $REPRESENTATIVES \
	--db $REPRESENTATIVES --self --id .84 \
	--iddef 1 --userout $OUTPUT \
	--userfields query+target+id \
	--maxaccepts 0 --query_cov .9 --maxhits 10

awk 'BEGIN{FS=";"}{print $1, $3, $5}' $OUTPUT > "${OUTPUT/.txt/_formatted.txt}"