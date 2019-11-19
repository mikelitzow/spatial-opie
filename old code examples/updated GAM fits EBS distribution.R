
library(mgcv)
library(dplyr)
library(ggplot2)

d = read.csv("ebs trawl taxa per Kotwicki & Lauth.csv")
d = d[which(d$bottom.temp > -10),]
d$logdepth = log(d$depth)

output = data.frame("year" = unique(d$year))
output_bin = data.frame("year" = unique(d$year))

spp.list = names(d)[2:23]

# covariates: lat, long, bottom.temp, depth, year, julian

for(i in 1:length(spp.list)) {
  
  # fit simple spatial GAM with field varying by year
  d$abund = d[,spp.list[i]]
  
  # fit delta-glm
  d$pres = ceiling(d$abund/1.0e10)
  model_bin = gam(pres ~ te(lat, long) + as.factor(year) +
                    bottom.temp + logdepth + julian + I(bottom.temp^2) +
                    I(logdepth^2) + I(julian^2), data=d, family="binomial")
  
  model_dens = gam(abund ~ te(lat, long) + as.factor(year) +
                     bottom.temp + logdepth + julian + I(bottom.temp^2) +
                     I(logdepth^2) + I(julian^2), data=d[which(d$abund>0),], family="Gamma")
  
  # only project to the same years in both datasets
  years = model_dens$xlevels$`as.factor(year)`[which(model_dens$xlevels$`as.factor(year)` %in% model_bin$xlevels$`as.factor(year)`)]
  
  df_new = expand.grid("lat" = seq(quantile(d$lat, 0.025),
                                   quantile(d$lat, 0.975), length.out=100),
                       "long" = mean(d$long),
                       "year" = as.numeric(years),
                       "bottom.temp" = mean(d$bottom.temp),
                       "logdepth" = mean(d$logdepth),
                       "julian" = mean(d$julian))
  
  # predict to new data frame
  df_new$predbin = predict(model_bin,
                           newdata = df_new, type="response")
  
  df_new$predpos = predict(model_dens,
                           newdata = df_new, type="response")
  df_new$pred = df_new$predbin * df_new$predpos # total
  
  # calculate weighted average of latitude, with weights = density
  g1 = group_by(df_new, year) %>%
    summarize(meanLat = weighted.mean(x = lat, w = pred),
              meanPres = mean(predbin))
  
  if(i==1) {
    output = data.frame("spp"= spp.list[i], "lat" = c(g1$meanLat), stringsAsFactors = FALSE)
    output_bin = data.frame("spp"= spp.list[i], "lat" = c(g1$meanPres), stringsAsFactors = FALSE)
  }
  if(i > 1) {
    output = rbind(output, data.frame("spp"= spp.list[i], "lat" = c(g1$meanLat), stringsAsFactors = FALSE))
    output_bin = rbind(output_bin, data.frame("spp"= spp.list[i], "lat" = c(g1$meanPres), stringsAsFactors = FALSE))
  }
  #g1 = ggplot(g1, aes(year, meanLat)) + geom_line() +
  #  geom_point() + ggtitle(spp.list[i])
  #print(g1)
}

library(ggplot2)
library(ggsidekick)
output = group_by(output, spp) %>%
  mutate(years=rev(rev(1982:2016)[1:length(lat)]))

ggplot(output, aes(x=years, y = lat)) + geom_line() + geom_point(size=0.4) +
  facet_wrap(~spp, scale="free_y") + theme_sleek() + xlab("Year") +
  ylab("Model - based center of gravity") + theme(strip.text.x = element_text(size = 4))

