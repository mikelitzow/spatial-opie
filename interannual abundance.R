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
dat$area <- as.factor(ifelse(DescTools::IsOdd(dat$CRUISE)==F, "NBS", "EBS"))

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
# structure...
# first, I'll confirm that check AICc values with and without the nugget in the correlation structure

model.comparison <- data.frame() # to catch the results

area <- as.factor(c("NBS", "EBS"))
group <- colnames(dat)[c(14:16,18,19)]

for(a in 1:length(area)) { # loop through each area

  # subset data for area of interest

  # THIS STEP DOESN'T WORK! NOT SURE WHY??
  # temp <- dat %>%
  #   filter(area==area[a]) 
  
  # use regular notation instead of dplyr...
  temp <- dat[dat$area==area[a],] 
  
  for(g in 1:length(group)) { # loop through each group

    # set formula for group of interest
    form <- as.formula(paste(group[g], "^0.25 ~ SURVEY_YEAR", sep=""))
    
    # record AICc values
    # m1 = with nugget, m2 = without nugget
    m1 <- gls(model=form, 
               correlation=corGaus(form=~LATITUDE + LONGITUDE | SURVEY_YEAR, 
                                   nugget = T), data=temp)
    m2 <- gls(model=form, 
              correlation=corGaus(form=~LATITUDE + LONGITUDE | SURVEY_YEAR, 
                                  nugget = F), data=temp)
    
    model.comparison <- rbind(model.comparison,
                              data.frame(area=area[a],
                                         group=group[g],
                                         nugget.T.AICc=AICc(m1),
                                         nugget.F.AICc=AICc(m2)))
    
  } # close group
} # close area

model.comparison # model with nugget is superior in every case

# so, re-run these loops and capture the annual abundance estimates for each year

abundance.estimates <- data.frame() # make a new object to catch results

for(a in 1:length(area)) { # loop through each area
  
  # subset data for area of interest
  temp <- dat[dat$area==area[a],] 
  
  for(g in 1:length(group)) { # loop through each group
    
    # set formula for group of interest
    form <- as.formula(paste(group[g], "^0.25 ~ SURVEY_YEAR", sep=""))
    
    # fit with nugget
    mod <- gls(model=form, 
              correlation=corGaus(form=~LATITUDE + LONGITUDE | SURVEY_YEAR, 
                                  nugget = T), data=temp)
   
    t <- summary(mod)$tTable
    
    abundance.estimates <- rbind(abundance.estimates,
                              data.frame(area=area[a],
                                         group=group[g],
                                         year=c(2010, 2017:2019),
                                         estimate=c(t[1,1], t[1,1]+t[2,1], t[1,1]+t[3,1], t[1,1]+t[4,1]),
                                         se=c(t[1,2], t[2,2], t[3,2], t[4,2]))) 
    
  } # close group
} # close area

abundance.estimates

# plot
theme_set(theme_bw())
# set a sex column
fem <- grep("FEM", abundance.estimates$group)
abundance.estimates$sex <- "Male"
abundance.estimates$sex[fem] <- "Female"
abundance.estimates$year <- as.factor(as.character(abundance.estimates$year))

ggplot(filter(abundance.estimates, sex=="Male"), aes(year, estimate)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymax=estimate+1.96*se, ymin=ifelse(estimate-1.96*se > 0, estimate-1.96*se, 0)), 
                position="dodge", width=0.5) +
  facet_grid(group~area, scales="free") +
  ylab("Fourth-root CPUE")

ggsave("figs/male annual abundance with spatial error.png", width=4, height=6, units='in')

# and females
ggplot(filter(abundance.estimates, sex=="Female"), aes(year, estimate)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymax=estimate+1.96*se, ymin=ifelse(estimate-1.96*se > 0, estimate-1.96*se, 0)), 
                position="dodge", width=0.5) +
  facet_grid(group~area, scales="free") +
  ylab("Fourth-root CPUE")

ggsave("figs/female annual abundance with spatial error.png", width=4, height=4, units='in')

# so, all of the interannual differences appear to change when errors are treated this way

###########################################
# for comparison's sake, look at lm results
abundance.estimates <- data.frame() # make a new object to catch results

for(a in 1:length(area)) { # loop through each area
  
  # subset data for area of interest
  temp <- dat[dat$area==area[a],] 
  
  for(g in 1:length(group)) { # loop through each group
    
    # set formula for group of interest
    form <- as.formula(paste(group[g], "^0.25 ~ SURVEY_YEAR", sep=""))
    
    # fit with nugget
    mod <- aov(formula=form, data=temp)

    t <- summary.lm(mod)$coefficients
    
    abundance.estimates <- rbind(abundance.estimates,
                                 data.frame(area=area[a],
                                            group=group[g],
                                            year=c(2010, 2017:2019),
                                            estimate=c(t[1,1], t[1,1]+t[2,1], t[1,1]+t[3,1], t[1,1]+t[4,1]),
                                            se=c(t[1,2], t[2,2], t[3,2], t[4,2]))) 
    
  } # close group
} # close area

abundance.estimates

# plot
theme_set(theme_bw())
# set a sex column
fem <- grep("FEM", abundance.estimates$group)
abundance.estimates$sex <- "Male"
abundance.estimates$sex[fem] <- "Female"
abundance.estimates$year <- as.factor(as.character(abundance.estimates$year))

ggplot(filter(abundance.estimates, sex=="Male"), aes(year, estimate)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymax=estimate+1.96*se, ymin=ifelse(estimate-1.96*se > 0, estimate-1.96*se, 0)), 
                position="dodge", width=0.5) +
  facet_grid(group~area, scales="free") +
  ylab("Fourth-root CPUE")

ggsave("figs/male annual abundance NO spatial error.png", width=4, height=6, units='in')

# and females
ggplot(filter(abundance.estimates, sex=="Female"), aes(year, estimate)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymax=estimate+1.96*se, ymin=ifelse(estimate-1.96*se > 0, estimate-1.96*se, 0)), 
                position="dodge", width=0.5) +
  facet_grid(group~area, scales="free") +
  ylab("Fourth-root CPUE")

ggsave("figs/female annual abundance NO spatial error.png", width=4, height=4, units='in')

#######################
# out of curiosity...compare with spatial models not including a nugget
abundance.estimates <- data.frame() # make a new object to catch results

for(a in 1:length(area)) { # loop through each area
  
  # subset data for area of interest
  temp <- dat[dat$area==area[a],] 
  
  for(g in 1:length(group)) { # loop through each group
    
    # set formula for group of interest
    form <- as.formula(paste(group[g], "^0.25 ~ SURVEY_YEAR", sep=""))
    
    # fit with nugget
    mod <- gls(model=form, 
               correlation=corGaus(form=~LATITUDE + LONGITUDE | SURVEY_YEAR, 
                                   nugget = F), data=temp)
    
    t <- summary(mod)$tTable
    
    abundance.estimates <- rbind(abundance.estimates,
                                 data.frame(area=area[a],
                                            group=group[g],
                                            year=c(2010, 2017:2019),
                                            estimate=c(t[1,1], t[1,1]+t[2,1], t[1,1]+t[3,1], t[1,1]+t[4,1]),
                                            se=c(t[1,2], t[2,2], t[3,2], t[4,2]))) 
    
  } # close group
} # close area

abundance.estimates

# plot
theme_set(theme_bw())
# set a sex column
fem <- grep("FEM", abundance.estimates$group)
abundance.estimates$sex <- "Male"
abundance.estimates$sex[fem] <- "Female"
abundance.estimates$year <- as.factor(as.character(abundance.estimates$year))

ggplot(filter(abundance.estimates, sex=="Male"), aes(year, estimate)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymax=estimate+1.96*se, ymin=ifelse(estimate-1.96*se > 0, estimate-1.96*se, 0)), 
                position="dodge", width=0.5) +
  facet_grid(group~area, scales="free") +
  ylab("Fourth-root CPUE")

ggsave("figs/male annual abundance with spatial error NO nugget.png", width=4, height=6, units='in')

# and females
ggplot(filter(abundance.estimates, sex=="Female"), aes(year, estimate)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymax=estimate+1.96*se, ymin=ifelse(estimate-1.96*se > 0, estimate-1.96*se, 0)), 
                position="dodge", width=0.5) +
  facet_grid(group~area, scales="free") +
  ylab("Fourth-root CPUE")

ggsave("figs/female annual abundance with spatial error NO nugget.png", width=4, height=4, units='in')

