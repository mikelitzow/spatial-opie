library(tidyverse)

# load data
dat <- read.csv("data/opie.data.csv")
range(dat$SURVEY_YEAR)

# check
head(dat)
str(dat)

unique(dat$SURVEY_YEAR)

# make a total CPUE colum
dat$total.cpue <- dat$TOTAL_FEMALE + dat$MALE_TOTAL

# look at the distribution
hist(dat$total.cpue^0.25) # fourth root -- not so good!
hist(log(dat$total.cpue+1))# log transformed - better

# for now, restrict to 1986-2019
dat <- dat %>%
  filter(SURVEY_YEAR >= 1986)

# identify stations sampled in every year
ff <- function(x) sum(!is.na(x))
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

# now...super-quickly
# will calculate average temperature in each year, 
# and average weighted by log CPUE

# first, change 0 total CPUE to NA
change <- dat$total.cpue==0
dat$total.cpue[change] <- NA

# and log-transform
dat$log.total.cpue <- log(dat$total.cpue)

raw.temp <- tapply(dat$NBT, dat$SURVEY_YEAR, mean, na.rm=T) # can impute NAs if needed!

ff <- function(x) weighted.mean(x, w=dat$log.total.cpue, na.rm=T)

occupied.temp <- NA
for(i in 1986:2019){
temp <- dat %>%
  filter(SURVEY_YEAR==i)

  occupied.temp[(i-1985)] <- weighted.mean(x=temp$NBT, w=temp$log.total.cpue, na.rm=T)
}

plot <- data.frame(year=1986:2019, raw.temp=raw.temp, occupied.temp=occupied.temp)

plot <- plot %>%
  gather(key, value, -year)

ggplot(plot, aes(year, value, color=key)) +
  theme_bw() +
  geom_line()

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
  geom_line()

cross.cor <- data.frame(year=1986:2019, diff=raw.temp-occupied.temp, raw.temp=raw.temp)

ccf(cross.cor$diff, cross.cor$raw.temp)

ggplot(cross.cor, aes(year, diff)) +
  theme_bw() +
  geom_line()

cross.cor <- cross.cor %>%
  gather(key, value, -year)
  
ggplot(cross.cor, aes(year, value, color=key)) +
  theme_bw() +
  geom_line()
