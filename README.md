# master-thesis

Supplementary material to the Master's thesis of Even Garvang. The scripts are for processing of metabarcoding data, based on the pipeline of [Frédéric Mahé](https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline)

A short description of each script follows below.

## Bioinformatics processing

### mergepairs.sh {- #mergepairs-sh}

Reads all the forward and reverse fastq-files in a directory and merges each paired-end sequence.

```bash
set -o errexit #exit on errors

# remove all loaded modules
module purge
# load the vsearch module
module load VSEARCH/2.9.1-foss-2018b

VSEARCH=$(which vsearch)
THREADS=16
ENCODING=33
JOB=merge_uncut

# list filenames containing R1

files=$(ls *-R1_001.fastq)

# loop over filenames and merge pairs

for file_R1 in $files; do
	
	# get filename for R2
	file_R2=${file_R1/R1/R2}
	samp="${file_R1/-R1*/}"
	
	#write info to logfile
	echo "Sample: ${samp}" >> ${JOB}.log

	# filter based on 18SV4 forward primer
  # Merge read pairs
  "${VSEARCH}" \
    --threads ${THREADS} \
    --fastq_mergepairs ${file_R1} \
    --reverse ${file_R2} \
    --fastq_ascii ${ENCODING} \
    --fastqout merged-fastq/${samp}-merged.fastq \
    --fastq_allowmergestagger \
    --quiet 2>> ${JOB}.log

done
```

### primer_demultiplex.sh {- #primer_demultiplex.sh}

Splits all files into two: one containing only sequences gained from the regular 18SV4 primers, and one containing those from the anti-metazoan primers.

```bash
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
```



### local_derep.sh {-}

Dereplicates each sample file, merging strictly identical sequences and annotating abundance.

```bash

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
```

### filter_negative.sh {-}

Reads the files of the negative controls, and removes any exact matches from the sample files.

```bash
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
```

### derep_cluster.sh {-}

Dereplicates all sequences globally, clusters with Swarm and checks for chimeras with VSEARCH's `--uchime_denovo` option.

```bash
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

```

### cut_pr2.sh {-}

Cuts the PR^2^ database to match our fragments with Cutadapt. This is an example with the general 18SV4 primers, a similar script was made with the anti-metazoan primers.

```bash
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
```

### filter_singletons.sh {-}

```bash
INPUT=common-cluster_1f_representatives.fas
OUTPUT=representatives_ns.fas

awk 'BEGIN{FS="[=;]"}NR%2==1&&$3>1{print;getline;print}' $INPUT > $OUTPUT

awk 'BEGIN{FS="[;]"}NR%2==1 {print substr($1, 2)}' \
  $OUTPUT > seeds_filtered

awk 'NR==FNR{id[$1]; next} $3 in id' seeds_filtered \
  ${INPUT/_representatives.fas/.stats} > ${OUTPUT/fas/stats}

awk 'BEGIN{FS=";"} NR==FNR{id[$1]; next} $1 in id' seeds_filtered \
  ${INPUT/_representatives.fas/.swarms} > ${OUTPUT/fas/swarms}

awk 'NR==FNR{id[$1]; next} $1 in id' seeds_filtered \
  ${INPUT/_representatives.fas/.struct} > ${OUTPUT/fas/struct}
```


### stampa.sh {-}

Queries the sequences against the cut PR^2^ database. Then it merges multiple best hits with the script stampa_merge.py and sorts by decreasing abundance. The script stampa_merge.py is written by Frédéric Mahé, available for download from: <https://github.com/frederic-mahe/stampa/raw/master/stampa_merge.py>

```bash
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
    --userout - | sed 's/;size=/_/ ; s/;//' > \
    common-tax/hits.representatives

module purge
# stampa_merge.py uses python 3
module load Python/3.8.6-GCCcore-10.2.0

# in case of multi-best hit, find the last-common ancestor
python scripts/stampa_merge.py $(pwd)/common-tax

# sort by decreasing abundance
sort -k2,2nr -k1,1d common-tax/results.representatives >\
  common-tax/common-ns-representatives.results
```


### OTU_contingency.sh {-}

Creates a contingency table of OTUs and samples. It calls the script OTU_contingency_table.py, which is written by Frédéric Mahé, available at <https://github.com/frederic-mahe/swarm/wiki/Fred's-metabarcoding-pipeline#build-the-otu-table>.

```bash
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
```

### self_match.sh {-}

Checks the pairwise similarities between all OTUs in  the data.

```bash
module purge
module load VSEARCH/2.9.1-foss-2018b

REPRESENTATIVES=common-swarms/representatives_ns.fas
OUTPUT=lulu/match_list_common.txt

vsearch --usearch_global $REPRESENTATIVES \
	--db $REPRESENTATIVES --self --id .84 \
	--iddef 1 --userout $OUTPUT \
	--userfields query+target+id \
	--maxaccepts 0 --query_cov .9 --maxhits 10

awk 'BEGIN{FS=";"}{print $1, $3, $5}' $OUTPUT >\
  "${OUTPUT/.txt/_formatted.txt}"
```

### lulu.R

Runs the LULU algorithm on an OTU table. Sums the two primers before curation, splits them afterwards.

```r
library(tidyverse)
library(lulu)
source("OTUtab-functions.R")

#### PRE-PROCESS OTU TABLE ####

OTU_raw <- read_OTU("data/common-OTU.table", split=TRUE)

OTU <- OTU_raw %>%
  filter(chimera == "N") %>%
  select(amplicon, contains("derep")) %>%
  column_to_rownames("amplicon")


cols <- str_match(colnames(OTU), "S[0-9][0-9]?") %>%
  unique()

# Sum the two samples from each primer, store the proportion between them

samp_prop <- matrix(nrow = nrow(OTU), ncol = ncol(OTU)/2, dimnames = list(rownames(OTU), cols))
x <- 1

for(j in seq(1, ncol(OTU), 2)){
  for(i in 1:nrow(OTU)){
    samp_prop[i,x] <- if(OTU[i,j+1]==0) 1 else OTU[i,j]/(OTU[i, j+1] + OTU[i,j])
  }
  x = x+1
}

OTU_sums = matrix(nrow = nrow(OTU), ncol = ncol(OTU)/2, dimnames = list(rownames(OTU), cols))
x <- 1

for(j in seq(1, ncol(OTU), 2)){
  for(i in 1:nrow(OTU)){
    OTU_sums[i,x] <- OTU[i,j] + OTU[i,j+1]
  }
  x = x+1
}

#### LULU-CURATION ####

matches <- read.table("data/match_list_common_formatted.txt")
curated <- lulu(as.data.frame(OTU_sums), matches)


lulu_table <- curated$curated_table %>%
  rownames_to_column("amplicon")

lulu_final <- OTU_raw %>%
  select(- contains("derep")) %>%
  right_join(lulu_table, by = "amplicon") %>%
  filter(chimera == "N")

write_tsv(lulu_final, "data/common2_lulu.tsv")
saveRDS(curated, "data/common2_lulu.rds")

curated <- readRDS("data/common2_lulu.rds")

#### SPLIT THE TWO PRIMER SETS ####

lulu_table_split <- matrix(
    nrow = nrow(lulu_table), 
    ncol = ncol(OTU),
    dimnames = list(rownames(curated$curated_table), colnames(OTU))) %>%
  as.data.frame()

samp_prop2 <- as.data.frame(samp_prop[rownames(curated$curated_table),])


x <- 1
for (j in seq(1, ncol(lulu_table_split), 2)){
  lulu_table_split[c(j, j+1)] <- cbind(curated$curated_table[x] * samp_prop2[x],
                                       curated$curated_table[x] * (1-samp_prop2[x]))
  x <- x+1
}

lulu_table_split2 <- rownames_to_column(round(lulu_table_split, 0), "amplicon")

lulu_final_split <- OTU_raw %>%
  select(- contains("derep")) %>%
  right_join(lulu_table_split2, by = "amplicon") %>%
  filter(chimera == "N")

write_tsv(lulu_final_split, "data/common2_lulu_split.tsv")
```

## Data analysis

### CTD-formatting.R {-}

Reads all CTD data from a folder and combines to one table. Requires that the file name structure is "YYYY-MM-DD.ID.cnv".

```r
# read in CTD data, write to file

library(oce)
library(ggplot2)
library(purrr)
library(stringr)
library(lubridate)

# create vector of filenames
ctd_files <- dir(path = "./CTD/", pattern = ".cnv$")

# read CTD data to a list of ctd class objects
ctd <- lapply(paste0("CTD/", ctd_files), read.ctd)

# function for converting ctd class objects to a data frame
# also extracting date+time, Brunt-vaisala, station ID from filenames
ctd_to_df <- function(ctd, filenames){
  ctd <- subset(ctd, sigmaT>10 & depth > 0)
  d <- as.data.frame(ctd@data)
  d$date_time <- ctd@metadata$startTime
  d$N2 <- swN2(ctd, derivs = "smoothing", df = 50)
  d$station_id <- str_match(filenames, "[.](.*)[.]")[,2]
  d$date <- as.Date(str_match(filenames, "^([^.]*)")[,1])

  return(d)
} 


# apply function to all objects and filenames
ctd_df <- purrr::map2(ctd, ctd_files, ctd_to_df)

# convert list of data frames to 1 data frame
ctd_full <- purrr::map_dfr(ctd_df, rbind)
```

### STD-formatting.R {-}

Manually reads .csv files and combines to a single table. Example shown for September data.

```r
library(tidyverse)
library(readxl)

samples <- read_xlsx("data/DNA_extraction_samples.xlsx") %>%
  mutate(Date = as.Date(date))

# function for reading STD file
# formats dates, creates depth from pressure, extracts
# only the relevant readings using the samples file.

read_STD <- function(filename, skip, filter_depth = 1){
  read_csv2(filename, skip = skip) %>%
    select(-c(X11, X12)) %>%
    mutate(
      Date = ifelse(str_length(Date) < 8,
                    paste0("0", Date),
                    Date), #dealing with annoying date shenanigans
      Date = as.Date(as.character(Date), format = "%d%m%Y"),
      Depth = marelac::sw_depth(Press/10, lat = 60)) %>%
    filter(Date %in% samples$Date,
           Press > filter_depth,
           Density > 10)
}

september <- read_STD("STD/2020-09.txt", skip = 17)

# infer station from depth
september %>% group_by(Ser) %>%
  summarize(max(Depth),
            max(Density))

september_id <- tribble(
  ~Ser, ~station_id,
  2, "DK1",
  3, "IM2",
  4, "FL1",
  5, "BF",
  6, "BN1"
)

# join station ID with the correct serials
september_full <- september %>%
  left_join(september_id)
```

### ctd-std-join.R

```r
library(tidyverse)

CTD <- read_tsv("CTD/CTD_full.tsv")
STD <- read_tsv("STD/STD_full.tsv")

# match names for easier combination
colnames(CTD)[c(4, 9, 10)] <- c("density", "oxygen_percent", "oxygen_mg_l")

# merge
std_ctd_merged <- dplyr::bind_rows(CTD, STD)

# write to file
write_tsv(std_ctd_merged, "STD_CTD_full.tsv")
```

### envvir-summary.R

```{r}
library(tidyverse)
library(lubridate)

std_ctd_merged <- read_tsv("STD_CTD_full.tsv")

#function for seconds since midnight
clockS = function(t){hour(t)*3600+minute(t)*60+second(t)}

ctd_summary <- std_ctd_merged %>%
  group_by(station_id, date) %>%
  summarise(max_depth = max(depth, na.rm = TRUE),
            pycnocline = max(unique(depth[N2 == max(N2, na.rm = TRUE)])),
            upr_temp = mean(temperature[depth<pycnocline], na.rm = TRUE),
            lwr_temp = mean(temperature[depth>pycnocline], na.rm = TRUE),
            upr_sal = mean(salinity[depth<pycnocline], na.rm = TRUE),
            lwr_sal = mean(salinity[depth>pycnocline], na.rm = TRUE),
            upr_dens = mean(density[depth<pycnocline], na.rm = TRUE),
            lwr_dens = mean(density[depth>pycnocline], na.rm = TRUE),
            time_start = min(clockS(date_time)))

write_tsv(ctd_summary, "ctd_envvar.tsv")
```

## Utility

Some utility functions are gathered in the script `OTUtab-functions.R`

```r

# contains some simple, useful functions for analyzing metabarcoding data
# Dependencies:
# tidyverse

library(tidyverse)
library(lubridate)

plot_reads <- function(data, group = sample, by_primer = TRUE){
  
  group_var <- enquo(group)
  
  baseplot <- ggplot(data, aes(!!group_var, reads)) + 
    geom_col() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  if(by_primer){
    return(baseplot + facet_wrap(~primer, nrow = 2))
  }
  return(baseplot)
  
}

plot_tax <- function(data, tax = division, by_primer = TRUE, position = "fill", xvar = sample){
  
  tax_var <- enquo(tax)
  x_var <- enquo(xvar)
  
  baseplot <- ggplot(data, aes(!!x_var, reads, 
                               fill = fct_lump(!!tax_var, n = 7))) + 
    geom_col(position = position) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  
  if(by_primer){
    return(baseplot + facet_wrap(~primer, nrow = 2))
  }
  return(baseplot)
}


count_OTU <- function(data, group, by_primer = TRUE){
  
  group_var <- enquo(group)
  
  if(by_primer){
    data %>% group_by(primer, !!group_var) %>%
      summarize(nOTU = n_distinct(OTU))
  }else{
    data %>% group_by(!!group_var) %>%
      summarize(nOTU = n_distinct(OTU))
  }
}


split_tax_PR2 <- function(data, 
                          tax_colname = taxonomy, 
                          tax = c("kingdom", "supergroup", "division", "class",
                                  "order", "family", "genus", "species")){
  
  tax_var <- enquo(tax_colname)
  
  data %>%
    separate(!!tax_var, into = tax, sep = "\\|")
}

lengthen_OTU <- function(data, pattern = "derep"){
  
  pivot_longer(data,
               cols = contains(pattern), 
               names_to = c("sample", "primer", NA),
               names_sep = "-",
               values_to = "reads")
}

read_OTU <- function(filename, split = FALSE){
  data <- read_tsv(filename) %>%
    mutate(OTU = paste0("OTU_", OTU)) %>%
    filter(chimera =="N")
  if(split) 
    return(split_tax_PR2(data))
  return(data)
}

write_fasta <- function(data, filename, OTU_col = "OTU", sequence_col = "sequence"){
  
  data <- as.data.frame(data)
  
  # overwrite if filename already exists
  cat("", file = filename, append = FALSE)
  
  # write to fasta
  cat(paste0(">OTU_", data[,OTU_col, drop = TRUE], "\n", data[,sequence_col, drop=TRUE], collapse = "\n"),
      file = filename,
      append = TRUE)
}

add_NA_date <- function(data, date, na_cols){
  new_data <- data
  new_data$date <- as.Date(date)
  new_data$month <- month(date)
  
  for(col in na_cols){
    new_data[col] <- NA
  }
  
  rbind(data, new_data)
}

#doesn't work for now
add_zerocount <- function(data, vars, zero_cols){
  
  data_full <- expand.grid(sapply(vars, function(col) unique(data[col])))
  names(data_full) <- vars
  d <- anti_join(data_full, data)
  for(col in zero_cols){
    d[col] <- 0
  }
  rbind(data, d)
}

geom_april <- function(log = FALSE){
  geom_polygon(data = data.frame(
    x=as.Date(c("2020-03-18","2020-03-18", "2020-04-27","2020-04-27")), 
    y = if(log) c(1e16,1e-16, 1e-16, 1e16) else c(Inf,-Inf, -Inf, Inf)),
    aes(x,y), inherit.aes = FALSE, fill = 1, alpha=0.2, )
}

geom_col_april <- function(log = FALSE){
  geom_polygon(data = data.frame(
    x=c(3.55,3.55, 4.45, 4.45), 
    y = if(log) c(1e16,1e-16, 1e-16, 1e16) else c(Inf,-Inf, -Inf, Inf)),
    aes(x,y), inherit.aes = FALSE, fill = 1, alpha=0.2)
}


manual_jitter <- function(coord, inc_scale = 0.5){
  if(length(coord) > 1){
    increment <- inc_scale*length(coord)/2
    new_coord <- seq(mean(coord) - increment, mean(coord) + increment, length.out = length(coord))
    return(new_coord)
  }
  return(coord)
}

plot_rar <- function(rar, xlab = "Sample size", ylab = "Species", 
                     title = "Rarefaction curve", col = alpha(rgb(0,0,0), 0.5), ...){
  
  by <- diff(attr(rar[[1]], "Subsample")[1:2])
  maxlen <- round(max(sapply(rar, function(x) max(length(x)))) * 1.03) * by
  maxspec <- round(max(sapply(rar, max)) * 1.05)
  
  
  plot(1:maxlen, seq(1, maxspec, length.out = maxlen), type = "n",
       xlab = xlab,
       ylab = ylab,
       las = 1)
  title(title, adj = 0)
  for(i in 1:length(rar)){
    endpos <- length(rar[[i]])
    lines(seq(1, endpos*by, by), rar[[i]], col = col, ...)
    graphics::text(endpos * by, rar[[i]][endpos], 
                   label = names(rar)[i],
                   pos = 4, offset = 0)
  }
}


calc_rar_slope <- function(rar, range = TRUE){
  
  by <- diff(attr(rar[[1]], "Subsample")[1:2])
  rar_slope <- rep(NA, length(rar))
  names(rar_slope) <- names(rar)
  
  for(i in 1:length(rar)){
    endpos <- length(rar[[i]])
    rar_slope[i] <- rar[[i]][(endpos - 100):(endpos-1)] %>%
      diff() %>%
      mean()
  }
  rar_slope <- rar_slope/by
  
  if(range) return(rar_slope[c(which.min(rar_slope), which.max(rar_slope))])
  
  return(rar_slope)
}


```