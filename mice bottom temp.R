# estimate missing values with mice!
require(mice)
require(Hmisc)
library(tidyverse)

# load data
dat <- read.csv("data/opie.data.csv")

# check
str(dat)

# restrict to stations sampled every year

# for now, restrict to 1988-2019 # changing this to 1988 start year - is that right?
dat <- dat %>%
  filter(SURVEY_YEAR >= 1988)

# make a total CPUE colum
dat$total.cpue <- dat$TOTAL_FEMALE + dat$MALE_TOTAL

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
sum(total==32) # 333! again, this has been changed to 1988-2019!

# restrict to these stations sampled each year
keepers <- total[total==32]
dat <- dat %>%
  filter(GIS_STATION %in% names(keepers))

# check for missing temperatures
temp <- tapply(dat$NBT, list(dat$SURVEY_YEAR, dat$GIS_STATION), mean)

# for now, I'll only include stations with no more than 22 NAs (at least 10 observations available)
j <- apply(temp, 2, function(x) sum(is.na(x)))
table(j)

keep <- apply(temp, 2, ff) >= 20

temp <- temp[,keep]

# Determine "best" predictor variables for each column
# first, get pairwise correlation coefficients
r <- rcorr(temp)$r # getting a NaN warning, but I don;t think this is a problem - don;t remember ever getting this before!

# Choose 20 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
pred <- t(apply(r,1, function(x) rank(1-abs(x))<=20)) # T for the 20 strongest correlations for each time series
diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE

Ni <- 100 # number of imputations
imp <- mice(temp, Ni, imp="norm", predictor = pred)

# Error in to - from : non-numeric argument to binary operator
#################################

#calculate mean and sd values from 100 imputations in final iteration
MI.sd <- MI.mean <- matrix(nrow = nrow(temp), ncol = ncol(temp))
dimnames(MI.sd) <- dimnames(MI.mean) <- dimnames(temp)

for(j in 1:ncol(temp)) { # go through each variable separately
  #j <- 2
  tt <- matrix(NA, nrow = nrow(temp), ncol = Ni)
  for(i in 1:Ni) { # and go through each imputation
    a <- complete(imp, i)
    tt[,i] <- a[,j]
  }
  MI.mean[,j] <- rowMeans(tt)
  MI.sd[,j] <- apply(tt, 1, sd)
}

# ok, now we need to replace the NAs!

mm <- is.na(d2$bottom.temp)
sum(mm) # 415 missing!

st.MI <- stack(as.data.frame(MI.mean))
st.MI$year <- 1982:2016
st.MI

d2 <- d2[order(d2$year),]
d2 <- d2[order(d2$station),]

identical(as.character(d2$station), as.character(st.MI$ind)) # TRUE!
identical(as.numeric(d2$year), as.numeric(st.MI$year)) # T!

d2$bottom.temp <- st.MI$values