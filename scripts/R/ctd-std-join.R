library(tidyverse)
library(lubridate)

CTD <- read_tsv("CTD/CTD_full.tsv")
STD <- read_tsv("STD/STD_full.tsv")

names(CTD)
names(STD)


# match names for easier combination

colnames(CTD)[c(4, 9, 10)] <- c("density", "oxygen_percent", "oxygen_mg_l")

std_ctd_merged <- dplyr::bind_rows(CTD, STD)


#write_tsv(std_ctd_merged, "STD_CTD_full.tsv")

# make giant plot

ggplot(std_ctd_merged, aes(-depth, density)) +
  geom_line(aes(col = station_id)) +
  facet_wrap(~month(date)) +
  coord_flip()

ggplot(std_ctd_merged, aes(-depth, N2)) +
  geom_line(aes(col = station_id)) +
  facet_wrap(~month(date)) +
  coord_flip()

# Get summary statistics for each station

#function for seconds since midnight
clockS = function(t){hour(t)*3600+minute(t)*60+second(t)}

ctd_summary <- std_ctd_merged %>%
  group_by(station_id, date) %>%
  summarise(max_depth = max(depth, na.rm = TRUE),
            pycnocline = max(unique(depth[N2 == max(N2, na.rm = TRUE)])),
            upr_temp = mean(temperature[depth<pycnocline], na.rm = TRUE),
            lwr_temp = mean(temperature[depth>pycnocline], na.rm = TRUE),
            upr10_temp = mean(temperature[depth<10], na.rm = TRUE),
            lwr10_temp = mean(temperature[depth>10], na.rm = TRUE),
            upr_sal = mean(salinity[depth<pycnocline], na.rm = TRUE),
            lwr_sal = mean(salinity[depth>pycnocline], na.rm = TRUE),
            upr10_sal = mean(salinity[depth<10], na.rm = TRUE),
            lwr10_sal = mean(salinity[depth>10], na.rm = TRUE),
            upr_dens = mean(density[depth<pycnocline], na.rm = TRUE),
            lwr_dens = mean(density[depth>pycnocline], na.rm = TRUE),
            time_start = min(clockS(date_time)))



sum(is.na(ctd_summary$pycnocline))

ggplot(ctd_summary, aes(pycnocline, station_id)) + geom_text(aes(label = factor(month(as.Date(date)))))

#write_tsv(ctd_summary, "ctd_envvar.tsv")


# new: add new data to summary
ctd_summary <- read_tsv("ctd_envvar.tsv")
niva <- readxl::read_xlsx("IO_2020_07_02_E09_binned.xlsx", sheet = "Data") 

niva <- niva %>%
  mutate(station_id = toupper(station_id)) %>%
  filter(station_id %in% ctd_summary$station_id | station_id == "EP1")

niva_summary <- niva %>% group_by(station_id) %>%
  summarize(upr10_temp = mean(Temperatur[Depth1<10], na.rm = TRUE),
            lwr10_temp = mean(Temperatur[Depth1>10], na.rm = TRUE),
            upr10_sal = mean(Saltholdighet[Depth1<10], na.rm = TRUE),
            lwr10_sal = mean(Saltholdighet[Depth1>10], na.rm = TRUE))

library(oce)
OF7 <- read.ctd("Sta0129.cnv")
plot(OF7)
OF7_df <- as_tibble(OF7@data)
names(OF7_df)[c(6, 7, 9)] <- c("salinity", "density", "what")
ggplot(OF7_df, aes(-pressure, temperature)) + geom_line() + coord_flip()

OF7_summary <- OF7_df %>%
  summarize(upr10_temp = mean(temperature[pressure<10], na.rm = TRUE),
            lwr10_temp = mean(temperature[pressure>10], na.rm = TRUE),
            upr10_sal = mean(salinity[pressure<10], na.rm = TRUE),
            lwr10_sal = mean(salinity[pressure>10], na.rm = TRUE),
            station_id = "IM2")

FL1_IM2_rows <- niva_summary[2,] %>%
  mutate(station_id = c("FL1"))

niva_new <- rbind(niva_summary, FL1_IM2_rows, OF7_summary) %>%
  mutate(date = as.Date("2020-07-02"))

ctd_summary_new <- bind_rows(ctd_summary, niva_new)

write_tsv(ctd_summary_new, "ctd_envvar_niva.tsv")

