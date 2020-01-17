# notes ----
# Time Series Figure
# Tyler Jackson
# Last Updated 2020-1-16

# load ----

library(tidyverse)
library(cowplot)
library(FNGr) # for tickr function from Ben Williams

# set x axis labels
x_axis <- tickr(data = tibble(yr = 1988:2019), yr, 3)

# data ----

# all time series
ts <- read_csv("./data/imputed_ts_bysex.csv")

# imputated temperature
imp_env <- read_csv("./data/imputed_temperature_timeseries.csv")

# figure ----

## create the plot of environmental timeseries
imp_env %>%
  ## mean bottom temperature
  select(1, 2) %>%
  ggplot(aes(x = SURVEY_YEAR, y = MEAN_BT, group = 1))+
  geom_point(color = "navyblue")+
  geom_line(color = "navyblue")+
  labs(y = "Mean Bottom \n Temperature (C)", x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(panel.grid = element_blank()) -> avg_bt

imp_env %>%
  ## cold pool areal extent
  select(1, 3) %>%
  ggplot(aes(x = SURVEY_YEAR, y = CPA, group = 1))+
  geom_point(color = "navyblue")+
  geom_line(color = "navyblue")+
  labs(y = bquote('Cold Pool Areal Extent ('~nm^2~')'), x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(panel.grid = element_blank()) -> cpa

imp_env %>%
  ## cold pool center of distribution
  select(1, 4) %>%
  ggplot(aes(x = SURVEY_YEAR, y = CP_COD, group = 1))+
  geom_point(color = "navyblue")+
  geom_line(color = "navyblue")+
  labs(y = "Cold Pool Center of \n Distribution (Latitude)", x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(panel.grid = element_blank()) -> cp_cod


## create the plot of biological timeseries
ts %>%
  # Temperature of Occupancy
  select(SURVEY_YEAR, SIZESEX, TEMP_OCC)  %>%
  ggplot(aes(x = SURVEY_YEAR, y = TEMP_OCC, group = SIZESEX, alpha = SIZESEX, shape = SIZESEX))+
  geom_point(color = "navyblue")+
  geom_line(color = "navyblue")+
  scale_alpha_manual(values = c(0.1, 0.1, 0.1, 0.1, 0.1, 1))+
  labs(y = "Temperature of \n Occupancy (C)", x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none") -> temp_occ

ts %>%
  # D95
  select(SURVEY_YEAR, SIZESEX, D95)  %>%
  ggplot(aes(x = SURVEY_YEAR, y = D95, group = SIZESEX, alpha = SIZESEX, shape = SIZESEX))+
  geom_point(color = "navyblue")+
  geom_line(color = "navyblue")+
  scale_alpha_manual(values = c(0.1, 0.1, 0.1, 0.1, 0.1, 1))+
  labs(y = bquote('D95 ('~nm^2~')'), x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none") -> d95

ts %>%
  # COD 
  select(SURVEY_YEAR, SIZESEX, LAT_COD) %>%
  ggplot(aes(x = SURVEY_YEAR, y = LAT_COD, group = SIZESEX, alpha = SIZESEX, shape = SIZESEX))+
  geom_point(color = "navyblue")+
  geom_line(color = "navyblue")+
  scale_alpha_manual(values = c(0.1, 0.1, 0.1, 0.1, 0.1, 1))+
  labs(y = "Crab Center of \n Distribution (Latitude)" , x = "")+
  scale_x_continuous(breaks = x_axis$breaks, labels = x_axis$labels)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = "none") -> cod


## combine plots
plot_grid(avg_bt, temp_occ, cpa, d95, cp_cod, cod, ncol = 2, align = "hv") -> ts_plot
 
## write plot
ggsave(filename = "./figs/all_timeseries.png", device = "png", width = 10, height = 8, 
       dpi = 300) 
