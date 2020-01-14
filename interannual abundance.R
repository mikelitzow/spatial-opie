library(tidyverse)
library(lme4)
library(nlme)
library(sp)
library(MuMIn)

# compare cpue for different sex/size categories in 2010, 2017-2019 in EBS & NBS

dat <- read.csv("data/opie.data.csv")
head(dat)

# limit to years of interest
dat <- dat %>%
  filter(SURVEY_YEAR %in% c(2010, 2017, 2018, 2019))

# separate into NBS/EBS
dat$area <- ifelse(DescTools::IsOdd(dat$CRUISE)==F, "NBS", "EBS")

# need to add 2018 NBS data, which are coded differently!
more.dat <- grep("NBS", dat$GIS_STATION)
dat$area[more.dat] <- "NBS"

# change year to factor!
dat$SURVEY_YEAR <- as.factor(as.character(dat$SURVEY_YEAR))

# first, look at NBS

NBS <- dat %>%
  filter(area=="NBS")

# look at the distributions... fourth-root transformation looks best
hist(NBS$MALE_31TO60)
hist(NBS$MALE_31TO60^0.25)
hist(log(NBS$MALE_31TO60+10, 10))


# Fit a model without correlated residuals
m0 <- gls(MALE_31TO60^0.25 ~ SURVEY_YEAR, data=NBS)
coordinates(NBS) <- c("LONGITUDE", "LATITUDE") 

# View semi variogram
for(i in c(2010, 2017, 2018, 2019)){
  select <- NBS$SURVEY_YEAR == i
  r <- resid(m0)[select]
  d <- as.vector(dist(coordinates(NBS)[select,]))
  p <- plot(Variogram(r, d), main = as.character(i))
  print(p)
} # Gaussian-ish?

m1 <- gls(MALE_31TO60^0.25 ~ SURVEY_YEAR, correlation=corGaus(form=~LATITUDE + LONGITUDE | SURVEY_YEAR, 
                                                                nugget = T), data=NBS)  
summary(m1)

m2 <- gls(MALE_31TO60^0.25 ~ SURVEY_YEAR, correlation=corGaus(form=~LATITUDE + LONGITUDE | SURVEY_YEAR, 
                                                                nugget = F), data=NBS) 
summary(m2)

AICc(m1, m2)

# are lat/long effects changing over years?
m3 <- gls(MALE_31TO60^0.25 ~ SURVEY_YEAR + LATITUDE:SURVEY_YEAR + LONGITUDE:SURVEY_YEAR,
            correlation=corGaus(form=~LATITUDE + LONGITUDE | SURVEY_YEAR, 
                                                                nugget = T), data=NBS)  
summary(m3)

m4 <- gls(MALE_31TO60^0.25 ~ SURVEY_YEAR + LATITUDE:SURVEY_YEAR + LONGITUDE:SURVEY_YEAR,
          correlation=corGaus(form=~LATITUDE + LONGITUDE | SURVEY_YEAR, 
                                                                nugget = F), data=NBS)

summary(m4)

AICc(m1, m2, m3, m4) 

# m3 is substantially the best
# but I wonder if it isn're more complicated than what we want here...don't know if we want to model spatial effects on
# abundance in each year?

# for now I'll just go through each area & size/sex group and model with the simple covariate 
# structure and the nugget in the correlation structure

areas <- c("NBS", "EBS")
group <- colnames(dat)[c(14:16,18,19)]

abundance.estimates <- data.frame() # make an object to catch results

for(a in 1:length(areas)) { # loop through each area
  
}