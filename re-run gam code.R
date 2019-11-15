library(mgcv)
library(dplyr)
library(ggplot2)

# re-run the script!
d2 <- read.csv("restricted stations EBS data.csv")

# drop the old rownumbers that were introduced into a column when the .csv file was created!
d2 <- d2[,-1]

# log transform depth
d2$logdepth = log(d2$depth)

output2 = data.frame("year" = unique(d2$year))
output_bin2 = data.frame("year" = unique(d2$year))

spp.list = names(d2)[2:23]

# covariates: lat, long, bottom.temp, depth, year, julian

for(i in 1:length(spp.list)) {
  
  # fit simple spatial GAM with field varying by year
  d2$abund = d2[,spp.list[i]]
  
  # fit delta-glm
  d2$pres = ceiling(d2$abund/1.0e10)
  model_bin = gam(pres ~ te(lat, long) + as.factor(year) +
                    bottom.temp + logdepth + julian + I(bottom.temp^2) +
                    I(logdepth^2) + I(julian^2), data=d2, family="binomial") 
  
  model_dens = gam(abund ~ te(lat, long) + as.factor(year) +
                     bottom.temp + logdepth + julian + I(bottom.temp^2) +
                     I(logdepth^2) + I(julian^2), data=d2[which(d2$abund>0),], family="Gamma")
  
  # only project to the same years in both datasets
  years = model_dens$xlevels$`as.factor(year)`[which(model_dens$xlevels$`as.factor(year)` %in% model_bin$xlevels$`as.factor(year)`)]
  
  df_new = expand.grid("lat" = seq(quantile(d2$lat, 0.025),
                                   quantile(d2$lat, 0.975), length.out=100),
                       "long" = mean(d2$long),
                       "year" = as.numeric(years),
                       "bottom.temp" = mean(d2$bottom.temp),
                       "logdepth" = mean(d2$logdepth),
                       "julian" = mean(d2$julian))
  
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
    output2 = data.frame("spp"= spp.list[i], "lat" = c(g1$meanLat), stringsAsFactors = FALSE)
    output_bin2 = data.frame("spp"= spp.list[i], "lat" = c(g1$meanPres), stringsAsFactors = FALSE)
  }
  if(i > 1) {
    output2 = rbind(output2, data.frame("spp"= spp.list[i], "lat" = c(g1$meanLat), stringsAsFactors = FALSE))
    output_bin2 = rbind(output_bin2, data.frame("spp"= spp.list[i], "lat" = c(g1$meanPres), stringsAsFactors = FALSE))
  }
  #g1 = ggplot(g1, aes(year, meanLat)) + geom_line() +
  #  geom_point() + ggtitle(spp.list[i])
  #print(g1)
}

library(ggplot2)
# library(ggsidekick)
output2 = group_by(output2, spp) %>%
  mutate(years=rev(rev(1982:2016)[1:length(lat)]))

# save output2
write.csv(output2, "GAM center of gravity restricted stations.csv")
# ggplot(output2, aes(x=years, y = lat)) + geom_line() + geom_point(size=0.4) +
#   facet_wrap(~spp, scale="free_y") + theme_sleek() + xlab("Year") +
#   ylab("Model - based center of gravity") + theme(strip.text.x = element_text(size = 4))


ggplot(output2, aes(x=years, y = lat)) + geom_line() + geom_point(size=0.4) +
  facet_wrap(~spp, scale="free_y") +  xlab("Year") +
  ylab("Model - based center of gravity") + theme(strip.text.x = element_text(size = 4))
