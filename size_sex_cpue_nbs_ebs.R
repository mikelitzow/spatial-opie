# notes ----
# Evaluate Abundance and CPUE by size sex group in EBS and NBS
# Tyler Jackson
# last updated: 2020/1/28

# load ----
library(tidyverse)
library(rsample)

## custom bootstrap function
f_boot_ci <- function(split){
  
  rsample::analysis(split) %>%
    mutate(abundance_mil = mean(cpue) * total_area / 1000000) %>%
    summarise(abundance_mil = mean(abundance_mil))
  
}

f_boot_ci2 <- function(split){
  
  rsample::analysis(split) %>%
    mutate(cpue = mean(cpue)) %>%
    summarise(cpue = mean(cpue))
  
}

# data ----

## NBS 
catch_nbs <- read_csv("./data/crab_nbs.csv")
haul_nbs <- read_csv("./data/haul_newtimeseries_nbs.csv")
strata_nbs <- read_csv("./data/strata_nbs_newtimeseries.csv")

## EBS
catch_ebs <- read_csv("./data/ebs_opilio_haul_data.csv")
haul_ebs <- read_csv("./data/haul_newtimeseries.csv")
strata_ebs <- read_csv("./data/ebs_opilio_strata.csv")


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
  mutate(size_sex = ifelse(sex == 1 & cw > 31 & cw <= 60, "male31to60",
                    ifelse(sex == 1 & cw > 61 & cw <= 90, "male61to90",
                    ifelse(sex == 1 & cw > 91 & cw <= 120, "male91to120",
                    ifelse(sex == 2 & clutch > 0, "mature_female",
                    ifelse(sex == 2 & clutch == 0, "immature_female", NA))))),
         area_swept = distance_fished * (net_width / 1000) * 0.539957^2) %>%
  group_by(year, Station, lat, lon, area_swept, size_sex) %>%
  summarise(num_crab = sum(sample_factor)) %>%
  right_join(strata, by = c("year", "Station")) %>%
  filter(!is.na(area_swept)) %>%
  pivot_wider(names_from = size_sex, values_from = num_crab) %>%
  pivot_longer(c(10:15), names_to = "size_sex", values_to = "num_crab") %>%
  filter(size_sex != "NA") %>%
  mutate(num_crab = replace_na(num_crab, 0),
         cpue = num_crab / area_swept,
         Region = ifelse(District == "NBS All", "NBS", "EBS")) %>%
  ungroup() -> cpue_long

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


# survey abundance estimates ----

## est abundance within region with naive error estimates by stratum_crab
## errors are summed among strata (assumes districts are independent)
cpue_long %>%
  group_by(year, Region, stratum_crab, size_sex) %>%
  summarise(district_area = mean(total_area),
            mean_cpue = mean(cpue),
            var_cpue = var(cpue),
            abundance_mil = district_area * mean_cpue / 1000000,
            var_abund_mil = district_area^2 * var_cpue / 1000000^2,
            stderr_mil = sqrt(var_abund_mil / n())) %>%
  group_by(year, Region, size_sex) %>%
  summarise(abundance_mil = sum(abundance_mil),
            stderr_mil = sum(stderr_mil)) %>%
  select(year, Region, size_sex, abundance_mil, stderr_mil) %>%
  mutate(lwr_95 = abundance_mil + qnorm(0.025) * stderr_mil,
         upp_95 = abundance_mil + qnorm(0.975) * stderr_mil) -> naive_est
### view plot
ggplot(naive_est, aes(x = factor(year), y = abundance_mil))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = lwr_95, ymax = upp_95), width = 0.3)+
  facet_wrap(~size_sex + Region, scales = "free_y", ncol = 2)

## est abundance within region with bootstrap error estimates by stratum_crab (same as District)
cpue_long %>%
  nest(-year, -District, -size_sex) %>%
  mutate(boot = map(data, ~rsample::bootstraps(., 1000))) %>% # bootstrap stations within strata
  unnest(boot) %>%
  mutate(models = map(splits, ~f_boot_ci(.x))) %>% # compute mean from every bootstrap sample
  unnest(models) %>%
  group_by(year, District, size_sex) %>%
  summarise(mean_abund = mean(abundance_mil),
            lwr_95 = quantile(abundance_mil, 0.025),
            upp_95 = quantile(abundance_mil, 0.975)) %>% # estimate 95% bootstrap CI with strata
  select(year, District, size_sex, mean_abund, lwr_95, upp_95) %>%
  mutate(Region = ifelse(District == "NBS All", "NBS", "EBS")) %>% # add region
  group_by(year, Region, size_sex) %>%
  summarise(mean_abund = sum(mean_abund),
            lwr_95 = sum(lwr_95),
            upp_95 = sum(upp_95)) %>% # sum abundances among regions
  full_join(naive_est %>%
              select(-lwr_95, -upp_95), 
            by = c("year", "Region", "size_sex")) %>% # replace bootstrap mean with point estimate
  select(year, Region, size_sex, abundance_mil, lwr_95, upp_95) -> boot_est
### view plot
ggplot(boot_est, aes(x = factor(year), y = abundance_mil))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = lwr_95, ymax = upp_95), width = 0.3)+
  facet_wrap(~size_sex + Region, scales = "free_y", ncol = 2)

## est abundance within region with weighted variance estimates by stratum_crab
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
         upp_95 = abundance_mil + qnorm(0.975) * stderr_mil) -> strat_est_crab
### view plot
ggplot(strat_est_crab, aes(x = factor(year), y = abundance_mil))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = lwr_95, ymax = upp_95), width = 0.3)+
  facet_wrap(~size_sex + Region, scales = "free_y", ncol = 2)

## est abundance within region with weighted variance estimates by stratum_crab
cpue_long %>%
  group_by(year, Region, stratum_fish, size_sex) %>%
  summarise(total_area = mean(total_area),
            mean_cpue = mean(cpue),
            var_cpue = var(cpue),
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
         upp_95 = abundance_mil + qnorm(0.975) * stderr_mil) -> strat_est_fish
### view plot
ggplot(strat_est_fish, aes(x = factor(year), y = abundance_mil))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = lwr_95, ymax = upp_95), width = 0.3)+
  facet_wrap(~size_sex + Region, scales = "free_y", ncol = 2)

# survey cpue estimates ----

## est cpue within region with naive variance estimates by stratum_crab
## errors are summed among strata (assumes districts are independent)
cpue_long %>%
  group_by(year, Region, stratum_crab, size_sex) %>%
  summarise(var_cpue = var(cpue),
            stderr = sqrt(var_cpue / n()),
            num_crab = sum(num_crab),
            area_swept = sum(area_swept)) %>%
  group_by(year, Region, size_sex) %>%
  summarise(cpue = sum(num_crab) / sum(area_swept),
            stderr = sum(stderr),
            lwr_95 = cpue + qnorm(0.025) * stderr,
            upp_95 = cpue + qnorm(0.975) * stderr) -> naive_est

## plot
ggplot(naive_est, aes(x = factor(year), y = cpue))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = lwr_95, ymax = upp_95), width = 0.3)+
  facet_wrap(~size_sex + Region, scales = "free_y", ncol = 2)

## est cpue within region with bootstrap error estimates by stratum_crab (same as District)
## does not perform well...
cpue_long %>%
  nest(-year, -District, -size_sex) %>%
  mutate(boot = map(data, ~rsample::bootstraps(., 1000))) %>% # bootstrap stations within strata
  unnest(boot) %>%
  mutate(models = map(splits, ~f_boot_ci2(.x))) %>% # compute mean from every bootstrap sample
  unnest(models) %>%
  group_by(year, District, size_sex) %>%
  summarise(mean_cpue = mean(cpue),
            lwr_95 = quantile(cpue, 0.025),
            upp_95 = quantile(cpue, 0.975)) %>% # estimate 95% bootstrap CI with strata
  select(year, District, size_sex, mean_cpue, lwr_95, upp_95) %>%
  mutate(Region = ifelse(District == "NBS All", "NBS", "EBS")) %>% # add region
  group_by(year, Region, size_sex) %>%
  summarise(mean_cpue = mean(mean_cpue),
            lwr_95 = mean(lwr_95),
            upp_95 = mean(upp_95)) %>% # sum abundances among regions
  full_join(naive_est %>%
              select(-lwr_95, -upp_95), 
            by = c("year", "Region", "size_sex")) %>% # replace bootstrap mean with point estimate
  select(year, Region, size_sex, cpue, lwr_95, upp_95) -> boot_est

## plot
ggplot(boot_est, aes(x = factor(year), y = cpue))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = lwr_95, ymax = upp_95), width = 0.3)+
  facet_wrap(~size_sex + Region, scales = "free_y", ncol = 2)

## est cpue within region with weighted variance estimates by stratum_crab
cpue_long %>%
  group_by(year, Region, stratum_crab, size_sex) %>%
  summarise(total_area = mean(total_area),
            num_crab = sum(num_crab),
            area_swept = sum(area_swept),
            var_cpue = var(cpue),
            stderr = sqrt(var_cpue / n())) %>%
  group_by(year, Region, size_sex) %>%
  mutate(bs_area = sum(total_area),
         w = total_area / bs_area) %>%
  summarise(cpue = sum(w * num_crab) / sum(w * area_swept),
            stderr = sum(w^2 * stderr),
            lwr_95 = cpue + qnorm(0.025) * stderr,
            upp_95 = cpue + qnorm(0.975) * stderr) -> strat_est_crab

## plot
ggplot(strat_est_crab, aes(x = factor(year), y = cpue))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = lwr_95, ymax = upp_95), width = 0.3)+
  facet_wrap(~size_sex + Region, scales = "free_y", ncol = 2)

## est cpue within region with weighted variance estimates by stratum_crab
cpue_long %>%
  group_by(year, Region, stratum_fish, size_sex) %>%
  summarise(total_area = mean(total_area),
            num_crab = sum(num_crab),
            area_swept = sum(area_swept),
            var_cpue = var(cpue),
            stderr = sqrt(var_cpue / n())) %>%
  group_by(year, Region, size_sex) %>%
  mutate(bs_area = sum(total_area),
         w = total_area / bs_area) %>%
  summarise(cpue = sum(w * num_crab) / sum(w * area_swept),
            stderr = sum(w^2 * stderr),
            lwr_95 = cpue + qnorm(0.025) * stderr,
            upp_95 = cpue + qnorm(0.975) * stderr) -> strat_est_fish

## plot
ggplot(strat_est_fish, aes(x = factor(year), y = cpue))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = lwr_95, ymax = upp_95), width = 0.3)+
  facet_wrap(~size_sex + Region, scales = "free_y", ncol = 2)