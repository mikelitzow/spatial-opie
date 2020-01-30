# notes ----
# Evaluate Abundance by  size sex group in EBS and NBS within 1mm bins
# Violin plot for NBS by year
# Tyler Jackson
# last updated: 2020/1/29

# load ----
library(tidyverse)
library(quantreg)

# data ----

## NBS 
catch_nbs <- read_csv("./data/Raw Survey Opilio Haul Data/nbs_opilio_haul_data.csv")
haul_nbs <- read_csv("./data/Raw Survey Opilio Haul Data/haul_newtimeseries_nbs.csv")
strata_nbs <- read_csv("./data/Raw Survey Opilio Haul Data/nbs_opilio_strata.csv")

## EBS
catch_ebs <- read_csv("./data/Raw Survey Opilio Haul Data/ebs_opilio_haul_data.csv")
haul_ebs <- read_csv("./data/Raw Survey Opilio Haul Data/haul_newtimeseries_ebs.csv")
strata_ebs <- read_csv("./data/Raw Survey Opilio Haul Data/ebs_opilio_strata.csv")


# data mgmt ----

## crab data NBS
catch_nbs %>%
  full_join(haul_nbs, by = c("HAULJOIN", "GIS_STATION")) %>%
  select(SURVEY_YEAR, CRUISE.x, GIS_STATION, WIDTH, SHELL_CONDITION, SEX, CLUTCH_SIZE,
         SAMPLING_FACTOR, DISTANCE_FISHED, NET_WIDTH, MID_LATITUDE, MID_LONGITUDE,
         HAUL_TYPE, PERFORMANCE) -> crab_nbs
names(crab_nbs) <- c("year", "cruise", "Station", "cw", "sc", "sex", "clutch", 
                     "sample_factor", "distance_fished", "net_width", "lat", "lon", 
                     "haul_type", "performance")

## crab data EBS
catch_ebs %>%
  full_join(haul_ebs, by = c("HAULJOIN", "GIS_STATION")) %>%
  filter(AKFIN_SURVEY_YEAR %in% c(2010, 2017:2019)) %>%
  select(AKFIN_SURVEY_YEAR, CRUISE.x, GIS_STATION, WIDTH, SHELL_CONDITION, SEX, CLUTCH_SIZE,
         SAMPLING_FACTOR, DISTANCE_FISHED.x, NET_WIDTH.x, MID_LATITUDE.x, MID_LONGITUDE.x,
         HAUL_TYPE.x, PERFORMANCE.x) -> crab_ebs
names(crab_ebs) <- c("year", "cruise", "Station", "cw", "sc", "sex", "clutch", 
                     "sample_factor", "distance_fished", "net_width", "lat", "lon", 
                     "haul_type", "performance")

## combine crab data
crab <- bind_rows(crab_nbs, crab_ebs)

## strata NBS
names(strata_nbs) <- c("Station", "District", "tows", "stratum", "total_area", "year")

## strata EBS
strata_ebs %>%
  filter(SURVEY_YEAR %in% c(2010, 2017:2019)) %>%
  select(STATION_ID, DISTRICT, TOWS, STRATUM, TOTAL_AREA_SQ_NM, SURVEY_YEAR) -> strata_ebs
names(strata_ebs) <- c("Station", "District", "tows", "stratum", "total_area", "year")

## combine strata
strata <- bind_rows(strata_nbs, strata_ebs)

## compute cpue by size-sex group for each station
crab %>%
  mutate(cw_bin = round(cw, 0),
         area_swept = distance_fished * (net_width / 1000) * 0.539957^2) %>%
  group_by(year, Station, lat, lon, area_swept, sex, cw_bin) %>%
  summarise(num_crab = sum(sample_factor)) %>%
  ungroup() %>%
  right_join(strata, by = c("year", "Station")) %>%
  filter(!is.na(area_swept)) %>%
  unite(col = size_sex, sex:cw_bin, sep = "_") %>%
  pivot_wider(names_from = size_sex, values_from = num_crab) %>%
  pivot_longer(grep("\\d", names(.)), names_to = "size_sex", values_to = "num_crab") %>%
  select(-"NA_NA") %>%
  mutate(num_crab = replace_na(num_crab, 0),
         cpue = num_crab / area_swept,
         Region = ifelse(District == "NBS All", "NBS", "EBS")) -> cpue_long

## Add groundfish stratum to cpue_long
## make 2018 EBS stratum_fish the same as 2018
haul_ebs %>%
  select(SURVEY_YEAR, GIS_STATION, STRATUM) %>%
  filter(SURVEY_YEAR %in% c(2010, 2017, 2019)) %>%
  bind_rows(haul_ebs %>%
              select(SURVEY_YEAR, GIS_STATION, STRATUM) %>%
              filter(SURVEY_YEAR %in% c(2019)) %>%
              mutate(SURVEY_YEAR = 2018)) %>%
  bind_rows(haul_nbs %>%
              select(SURVEY_YEAR, GIS_STATION, STRATUM)) %>%
  rename(year = SURVEY_YEAR,
         Station = GIS_STATION,
         stratum_fish = STRATUM) %>%
  right_join(cpue_long, by = c("year", "Station")) %>%
  rename(stratum_crab = stratum) %>%
  mutate(stratum_crab = ifelse(stratum_crab == 999, 10, stratum_crab)) -> cpue_long

# stratified estimate of abundance in 1 mm bins ----

cpue_long %>%
  group_by(year, Region, stratum_crab, size_sex) %>%
  summarise(total_area = mean(total_area),
            mean_cpue = mean(cpue),
            var_cpue = var(cpue),
            stratum_stations = n(),
            abundance_mil = mean(total_area) * mean_cpue / 1000000,
            var_abund_mil = mean(total_area)^2 * var_cpue / 1000000^2,
            stderr_mil = sqrt(var_abund_mil / n())) %>%
  group_by(year, Region, size_sex) %>%
  mutate(bs_area = sum(total_area),
         w = total_area / bs_area) %>% # define weights
  summarise(abundance_mil = sum(w * abundance_mil),
            stderr_mil = sum(w^2 * stderr_mil)) %>%
  select(year, Region, size_sex, abundance_mil, stderr_mil) %>%
  mutate(lwr_95 = abundance_mil + qnorm(0.025) * stderr_mil,
         upp_95 = abundance_mil + qnorm(0.975) * stderr_mil) -> stratum_est_crab

# violin plot of size comp ----

stratum_est_crab %>%
  separate(size_sex, into = c("sex", "size"), sep = "_") %>%
  mutate(size = as.numeric(size),
         sex = ifelse(sex == "1", "Male", "Female"),
         sex = factor(sex, levels = c("Male", "Female"))) %>%
  filter(size >= 30, 
         abundance_mil != 0,
         Region == "NBS") %>%
  group_by(year, sex) %>%
  mutate(w = abundance_mil / sum(abundance_mil),
         median = spatstat::weighted.median(size, w = w),
         q90 = spatstat::weighted.quantile(size, w = w, probs = 0.90)) %>%
  ggplot()+
  geom_violin(aes(x = factor(year), y = size, weight = w), adjust = 0.3, fill = "grey80")+
  geom_boxplot(aes(x = factor(year), y = size, weight = w), width = 0.1, outlier.alpha = 0.5)+
  #geom_point(aes(x = factor(year), y = median), size = 3)+
  #geom_point(aes(x = factor(year), y = q90), shape = 4, size = 3)+
  labs(x = NULL, y = "Carapce Width (mm)")+
  scale_fill_manual(values = c("grey40"))+
  scale_y_continuous(breaks = seq(30, 150, 15))+
  facet_wrap(~sex, labeller = )+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggsave("file_name.png", width = 4, height = 4, units = "in")
