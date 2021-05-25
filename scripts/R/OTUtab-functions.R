# contains some simple, useful functions for analyzing metabarcoding data

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
