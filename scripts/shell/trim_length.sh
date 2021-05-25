module purge
module load cutadapt/2.7-GCCcore-8.3.0-Python-3.7.4

for file in *.fastq; do
	cutadapt --minimum-length 300 \
	--maximum-length 450 $file > common-trim/$file 2>> common-trim/trim.log

done