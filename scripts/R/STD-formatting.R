library(tidyverse)
library(readxl)
library(marelac)

samples <- read_xlsx("data/DNA_extraction_samples.xlsx") %>%
  mutate(Date = as.Date(date))

read_STD <- function(filename, skip, filter_depth = 1){
  read_csv2(filename, skip = skip) %>%
    select(-c(X11, X12)) %>%
    mutate(Date = ifelse(str_length(Date) < 8, paste0("0", Date), Date), #dealing with annoying date shenanigans
           Date = as.Date(as.character(Date), format = "%d%m%Y"),
           Depth = marelac::sw_depth(Press/10, lat = 60)) %>%
    filter(Date %in% samples$Date,
           Press > filter_depth,
           Density > 10)
}

# start with STD
august <- read_STD("STD/2020-08.txt", skip = 16) %>%
  filter(Ser != 4)

august %>%
  group_by(Ser) %>%
  summarize(max(Depth),
            max(Density))

# Getting station info from station order
# although something seems to be wrong with depth of FL1?

august_id <- tribble(
  ~Ser, ~station_id,
  3, "BN1",
  5, "BF",
  6, "DK1",
  7, "FL1",
  8, "IM2"
)

august_full <- august %>%
  left_join(august_id)

ggplot(august_full, aes(-Press, Density)) +
  geom_line() +
  coord_flip() +
  facet_wrap(~station_id)

# SEPTEMBER

september <- read_STD("STD/2020-09.txt", skip = 17)

september %>% group_by(Ser) %>%
  summarize(max(Depth),
            max(Density))

# infer station from depth

september_id <- tribble(
  ~Ser, ~station_id,
  2, "DK1",
  3, "IM2",
  4, "FL1",
  5, "BF",
  6, "BN1"
)

september_full <- september %>%
  left_join(september_id)

ggplot(september_full, aes(-Press, Density)) +
  geom_line() +
  coord_flip() +
  facet_wrap(~station_id)

# OCTOBER

october <- read_STD("STD/2020-10.txt", skip = 18)

october %>%
  group_by(Ser) %>%
  summarize(max(Depth),
            max(Density))

#infer from max depth

october_id <- tribble(
  ~Ser, ~station_id,
  10, "BF",
  11, "BN1",
  12, "DK1",
  13, "FL1",
  14, "IM2"
)

october_full <- left_join(october, october_id)

ggplot(october_full, aes(-Press, Density)) +
  geom_line() +
  coord_flip() +
  facet_wrap(~station_id)


# Nov-dec

nov_dec <- read_STD("STD/2020-11.12.txt", skip = 29)

nov_dec %>% group_by(Ser) %>%
  summarize(max(Depth),
            max(Density))


nov_dec_id <- tribble(
  ~Ser, ~station_id,
  3, "BN1",
  4, "BF",
  5, "DK1",
  6, "FL1",
  7, "IM2",
  18, "BN1",
  19, "BF",
  20, "DK1",
  21, "IM2",
  22, "FL1"
)

nov_dec_full <- left_join(nov_dec, nov_dec_id)




STD_full <- rbind(
  august_full,
  september_full,
  october_full,
  nov_dec_full
)

col_names <- c("serial_no", "measurement_no", "salinity", "temperature",
               "oxygen_percent", "oxygen_mg_l", "density", "pressure",
               "date", "time", "depth", "station_id")


colnames(STD_full) <- col_names

std_list <- split(STD_full, list(STD_full$date, STD_full$station_id))
std_list_N2 <- lapply(std_list, function(x){ 
  x$N2 <- oce::swN2(x$pressure, x$density, df = 50)
  return(x)
  })

std_final <- purrr::map_dfr(std_list_N2, rbind) %>%
  mutate(date_time = as.POSIXct(paste(date, time)))


write_tsv(std_final, "STD/STD_full.tsv")
