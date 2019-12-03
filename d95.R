library(tidyverse)
library(maps)
library(mapdata)

# try to code up a selection of the stations containing 95% of CPUE in each year of the survey

# load data
dat <- read.csv("data/opie.data.csv")
range(dat$SURVEY_YEAR)

# check
head(dat)
str(dat)

unique(dat$SURVEY_YEAR)

# make a total CPUE colum
dat$total.cpue <- dat$TOTAL_FEMALE + dat$MALE_TOTAL

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

# now...loop through and select the stations making up 95% of cpue in a year
years <- unique(dat$SURVEY_YEAR) # get an object of all the years in the data
years <- sort(years) # put them in order to keep things clean

output <- data.frame() # make an object to hold results

for(i in 1:length(years)){ # loop through each year
 # i <- 1
  temp <- dat %>%
    filter(SURVEY_YEAR==years[i]) # select the data 
  
  temp <- temp %>%
    arrange(desc(total.cpue)) # sort by cpue (large:small)
  
  temp$proportion.cpue <- temp$total.cpue / sum(temp$total.cpue) # calculate the proportion of total cpue for each station
  
  temp$cumulative.cpue <- cumsum(temp$proportion.cpue) # get the cumulative proportion for each station
  
  # and save results!
  output <- rbind(output, 
                  data.frame(year=years[i], 
                             station=temp$GIS_STATION,
                             d95 = temp$cumulative.cpue <= 0.95))
}
