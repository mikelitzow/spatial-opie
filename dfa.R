library(tidyverse)
library(MARSS)
library(broom)

dat <- read.csv("data/annual.environmental.data.csv", row.names = 1)

# select the columns we want, and then transpose and turn into a matrix 
# (these are requirements for MARSS - column names need to be time steps)

dat <- dat %>%
  select(AVG_BT, CP_EXTENT, CP_COD)

# look at the distributions quickly!
plot.dat <- dat %>%
  gather()

theme_set(theme_bw())

ggplot(plot.dat, aes(value)) +
  geom_histogram() +
  facet_wrap(~key, scales="free")
  
# not terrible?

# now set up for MARSS
dat <- as.matrix(t(dat))

# too easy!

# set up forms of R matrices - these are the four candidate error structures
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")

# I'm using the default convergence criteria!

# make an object to save output
model.data = data.frame()

# fit models & store results
for(R in levels.R) {
  for(m in 1) {  # allowing only one shared trend as we only have 3 time series!
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(dat, model=dfa.model,
                 form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

# calculate delta-AICc scores, sort in descending order, and compare
model.data$dAICc <- model.data$AICc-min(model.data$AICc)
model.data <- model.data %>%
  arrange(dAICc)
model.data

# unconstrained errors the best!

# save the model comparison table
dir.create("output", showWarnings = FALSE)

write.csv(model.data, "output/environment dfa model selection.csv")

# now run the best model, save the trend estimates and CIs, and plot loadings/trend
model.list = list(A="zero", m=1, R="unconstrained")
mod <- MARSS(dat, model=model.list, z.score=TRUE, form="dfa")


# get CI and plot loadings...
modCI <- MARSSparamCIs(mod)
modCI

plot.CI <- data.frame(names=rownames(dat), mean=modCI$par$Z, upCI=modCI$par.upCI$Z,
                           lowCI=modCI$par.lowCI$Z)


plot.CI$names <- reorder(plot.CI$names, plot.CI$mean)

# don't know if I really need this step - leftover from previous interations!
dodge <- position_dodge(width=0.9)

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

dir.create("figs", showWarnings = FALSE)

ggplot(plot.CI, aes(x=names, y=mean)) + 
  geom_bar(position="dodge", stat="identity", fill=cb[2]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1,  size=12), axis.title.x = element_blank()) +
  geom_hline(yintercept = 0) + 
  ggtitle("DFA loadings")

ggsave("figs/environment dfa loadings.png", width=3, height=6, units='in')

# and the shared trend!
trend.CI <- tidy(mod, type="states")
trend.CI$year <- 1988:2019

ggplot(trend.CI, aes(x=year, y=estimate)) + 
  geom_line(color=cb[3]) +
  geom_ribbon(aes(ymax=conf.high, ymin=conf.low), fill=cb[3], alpha=0.2) +
  ylab("Shared trend") + 
  theme(axis.title.x = element_blank()) +
  geom_hline(yintercept = 0) + 
  ggtitle("DFA trend")

ggsave("figs/environment dfa trend.png", width=6, height=5, units='in')

# and save the trend output
write.csv(trend.CI, "output/environment dfa trend.csv", row.names = F)
