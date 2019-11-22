library(tidyverse)
library(maps)
library(mapdata)

# initial / exploratory analysis for analysis of changing opilio distributions
# in the EBS / NBS

# two general questions that I want to address:
# 1) how have avg. bottom temperatures changed in the EBS?
# 2) how have bottom temperatures where oplio occur in the EBS changed?

# the idea I'm interested in here is whether opies have successfully buffered 
# climate change effects by moving as temperature changes, or whether the changes in 
# temp have been too rapid / extreme to keep up with

# load data
dat <- read.csv("data/opie.data.csv")
range(dat$SURVEY_YEAR)

# check
head(dat)
str(dat)

unique(dat$SURVEY_YEAR)

# make a total CPUE colum
dat$total.cpue <- dat$TOTAL_FEMALE + dat$MALE_TOTAL

# look at the CPUE distribution
hist(dat$total.cpue^0.25) # fourth root -- not so good!
hist(log(dat$total.cpue+1))# log transformed - better

# for now, restrict to 1986-2019
dat <- dat %>%
  filter(SURVEY_YEAR >= 1986)

# identify stations sampled in every year
ff <- function(x) sum(!is.na(x)) # create a function that counts the non-missing cases
samples <- tapply(dat$total.cpue, list(dat$SURVEY_YEAR, dat$GIS_STATION), ff)
samples

# change NA to 0
change <- is.na(samples)
samples[change] <- 0

# and look at total # of samples
total <- colSums(samples)
total

# how many stations are sampled each year?
sum(total==34) # 305!

# restrict to these stations sampled each year
keepers <- total[total==34]
dat <- dat %>%
  filter(GIS_STATION %in% names(keepers))

# and check 
check <- data.frame(one=unique(dat$GIS_STATION)[order(unique(dat$GIS_STATION))], two=names(keepers))
check # looks good!

# now plot the mean lat/long for these 305 stations sampled in every year
mean.lat <- tapply(dat$LATITUDE, as.character(dat$GIS_STATION), mean) 
# changing lat to a character allows us to drop missing values that are included when lat is a factor!
mean.long <- tapply(dat$LONGITUDE, as.character(dat$GIS_STATION), mean) 

plot(360+mean.long, mean.lat) # change longitude to ºE to plot
map('world2Hires', 'usa', fill=T, add=T, 
    lwd=0.5, col="lightyellow3")


# now...super-quickly
# will calculate average temperature in each year, 
# and average weighted by log CPUE

# first, change 0 total CPUE to NA
change <- dat$total.cpue==0
dat$total.cpue[change] <- NA

# and log-transform
dat$log.total.cpue <- log(dat$total.cpue)

# check for NA temps
ff <- function(x) sum(is.na(x))
temp.na <- tapply(dat$NBT, list(dat$SURVEY_YEAR, as.character(dat$GIS_STATION)), ff)
View(temp.na) # tons of NAs!

raw.temp <- tapply(dat$NBT, dat$SURVEY_YEAR, mean, na.rm=T) # can impute NAs if needed!

ff <- function(x) weighted.mean(x, w=dat$log.total.cpue, na.rm=T) # new function to calculate mean temp weighted by log cpue

occupied.temp <- NA
for(i in 1986:2019){
temp <- dat %>%
  filter(SURVEY_YEAR==i)

  occupied.temp[(i-1985)] <- weighted.mean(x=temp$NBT, w=temp$log.total.cpue, na.rm=T)
}

plot <- data.frame(year=1986:2019, raw.temp=raw.temp, occupied.temp=occupied.temp) # put raw and occupied temp into data frame to plot

plot <- plot %>%
  gather(key, value, -year)

ggplot(plot, aes(year, value, color=key)) +
  theme_bw() +
  geom_line()

# very little difference between the two time series! 
# try with non-transformed cpue

occupied.temp <- NA
for(i in 1986:2019){
  temp <- dat %>%
    filter(SURVEY_YEAR==i)
  
  occupied.temp[(i-1985)] <- weighted.mean(x=temp$NBT, w=temp$total.cpue, na.rm=T)
}

plot <- data.frame(year=1986:2019, raw.temp=raw.temp, occupied.temp=occupied.temp)

plot <- plot %>%
  gather(key, value, -year)

ggplot(plot, aes(year, value, color=key)) +
  theme_bw() +
  geom_line() +
  ylab("ºC")

# so there seems that there might be buffering during 2014-2017, less so in 2018-2019

ggsave("plots/raw and occupied temp.png", width = 6, height = 4, units='in')

# check for lags
cross.cor <- data.frame(year=1986:2019, diff=raw.temp-occupied.temp, raw.temp=raw.temp, occupied=occupied.temp)

ccf(cross.cor$diff, cross.cor$raw.temp) # no meaningful correlations at various lags
ccf(cross.cor$occupied, cross.cor$raw.temp) # strong correlation at lag 0

# plot the difference time series
ggplot(cross.cor, aes(year, diff)) +
  theme_bw() +
  geom_line()

cross.cor <- cross.cor %>%
  select(-occupied) %>%
  gather(key, value, -year)
  
ggplot(cross.cor, aes(year, value, color=key)) +
  theme_bw() +
  geom_line() +
  ylab("ºC")
ggsave("plots/raw temp and diff.png", width = 6, height = 4, units='in')
