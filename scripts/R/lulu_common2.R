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