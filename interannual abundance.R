library(tidyverse)
library(lme4)
library(nlme)

# compare cpue for different sex/size categories in 2010, 2017-2019 in EBS & NBS

dat <- read.csv("data/opie.data.csv")
head(dat)

# create NBS/EBS code!

dat$area <- ifelse(DescTools::IsOdd(dat$CRUISE)==T, "EBS", "NBS")

# limit to years of interest
dat <- dat %>%
  filter(SURVEY_YEAR %in% c(2010, 2017:2019))

# change year to factor!
dat$SURVEY_YEAR <- as.factor(as.character(dat$SURVEY_YEAR))

# first, look at NBS

NBS <- dat %>%
  filter(area=="NBS")

hist(NBS$MALE_31TO60)

# better transform!

m1 <- gls(MALE_31TO60^0.25 ~ SURVEY_YEAR, correlation=corGaus(form=~LATITUDE / LONGITUDE), data=NBS)
summary(m1)
