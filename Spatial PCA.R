# Spatial and temporal patterns of total CPUE and CPUE of 
# selected taxa of groundfish in the Eastern Bering Sea

setwd("C:/Dropbox/Research/Projects/BSIERP/Analysis/Fish/Spatial patterns CPUE")

load("Spatial.RData")

require(mgcv)
require(lattice)
require(maps)
require(mapdata)

#################################################################################
# Multivariate analyses of spatial patterns in species CPUE
# First, select stations with no missing values and save average lat/long
# (exact locations differ among years):

# Also check for outliers (unusual lat/long in a given year):
lat <- tapply(BShauls$LATITUDE, list(BShauls$YEAR, BShauls$STATION), mean, na.rm=T)
long <- tapply(BShauls$LONGITUDE, list(BShauls$YEAR, BShauls$STATION), mean, na.rm=T)
j <- apply(lat, 2, function(x) sum(is.na(x)))
table(j)  # Number of stations with n missing values (years)

 # Selects only stations that were always sampled!
lat <- lat[,j==0]   
long <- long[,j==0]
st <- dimnames(lat)[[2]]  # SAVE stations!
yrs <- dimnames(lat)[[1]]
identical(st, dimnames(long)[[2]])  # check

# Largest absolute deviations in latitude and longitude from the mean position of a given station
hist(apply(lat, 2, function(x) max(x - mean(x))))
hist(apply(long, 2, function(x) max(x - mean(x))))

# Save mean lat/long:
lat <- apply(lat,2,mean)
long <- apply(long, 2, mean) 



####################################################################################################
# Patterns by species:

# Set up matrices and vectors for output for various spatial measures:
PC1 <- CG.lat <- CG.long <- Inertia <- avg.T <- div.CPUE <- matrix(NA, length(yrs), ncol(BScpue))
n <- dimnames(BScpue)[[2]]
dimnames(PC1) <- dimnames(avg.T) <- dimnames(div.CPUE) <- list(yrs, n)
dimnames(CG.lat) <- dimnames(CG.long) <- dimnames(Inertia) <- list(yrs, n)
EigVec <-  matrix(NA, length(st), ncol(BScpue))
dimnames(EigVec) <- list(st, n)
Var.expl <- matrix(NA, ncol(PC1), 5); dimnames(Var.expl) <- list(n, paste("PC", 1:5, sep=""))
Dev.expl <- numeric(ncol(PC1)); names(Dev.expl) <- n



#######################################################################
# Spatial PCA of species cpue by station

# The code below uses a number of functions and objects from othe sources:
# several function that are defined in a separate script (from M. Woillez)
# contours for the eastern Bering Sea!
# Object 'cpue.ann' - Need to run 
source("Spatial_indicators_functions.r")  # Make sure the file is in current working directory!
load("EBS.RData")
map.EBS()  # check to see if it works, should get "pretty" map
source("CPUE - Stratified random estimate.R")

# Graphical output:
win.graph(9,7.5)

# To save graphical output to pdf file, open pdf file first (close after running for loop using 'dev.off()'):
pdf("SpatialPCA.pdf", width=9, height=7.5)

par(omi=c(0,0,0.5,0), mar=c(4.2,4.2,2,1))
layout(matrix(c(2,3,0,1),2,2), width = c(2,1), height=c(1.5,1))
# layout.show(3)

years <- as.numeric(yrs)        


# Highlight and run entire for loop:

for(i in 1:ncol(BScpue)) {
 print(Sys.time())   # keep track of how long each loop takes!
# Select CPUE for given species and tabulate mean by year/station
 sp <- dimnames(BScpue)[[2]][i]
 cpue <- BScpue[,sp]
 cpue <- tapply(cpue, list(BShauls$YEAR, BShauls$STATION), mean, na.rm=T)
 cpue <- cpue[, st]  # Select stations sampled each year
 if(sum(is.na(cpue)) > 0) stop("Missing station/year combination in cpue")

# Eliminate stations where species was never caught:
 j <- apply(cpue,2,sum) > 0
 cpue <- cpue[,j]
 sp.lat <- lat[j]
 sp.long <- long[j]

# Average temperature across area of distribution:
 j <- is.element(BShauls$STATION, dimnames(cpue)[[2]])
 avg.T[,i] <- tapply(BShauls$GEAR_TEMP[j], BShauls$YEAR[j], mean, na.rm=T)

# Compute and save diversity of catches across stations: 
 f <- function(x) {
  x <- x[x > 0]
  pi <- x/sum(x)
  - sum(pi * logb(pi, base = 2))
 }
 div.CPUE[,i] <- apply(cpue,1,f) 

# Compute center of gravity and Inertia (Woillez et al 2009)
# See 'Spatial_indicators_functions.R'
 for(j in 1:nrow(cpue)) {
    CGI <- cgi(sp.long, sp.lat, cpue[j,], plot=F)
    CG.long[j,i] <- CGI$xcg
    CG.lat[j,i] <- CGI$ycg
    Inertia[j,i] <- CGI$I
 }

# Remove annual means across stations (CPUE anomalies)
 anom <- sweep(cpue^0.25, 1, apply(cpue^0.25,1,mean))
# Singular Value Decomposition
 anom.svd <- svd(var(anom))
# Extract and save first eigenvector, change sign for consistency to get
# increasing trend over time:
 evec <- anom.svd$u[,1]
# Time series of first PC:
 pc1 <- as.vector(scale(anom %*% evec))
 cf <- coef(lm(pc1 ~ years))[2]
 if(cf < 0) evec <- -evec
 EigVec[dimnames(cpue)[[2]],i] <- evec
# (Final) time series of first PC:
 pc1 <- as.vector(scale(anom %*% evec))

# Extract and plot eigenvalues:
 ev <- anom.svd$d
 x <- barplot(ev[1:8]/sum(ev)*100, xlab = "Principal Component", ylab="Proportion of variance")
 axis(1, at=x, 1:8, col=0)
 title("Scree plot")
# Save variance explained by first 5 eigenvectors
 Var.expl[i,] <- ev[1:5]/sum(ev) * 100  # Variance explained (%)

# Spatial pattern:
 K <- min(trunc(0.6*ncol(anom)), 60) 
 fit <- gam(evec ~ s(sp.lat, sp.long, k=K))
 summary(fit)
 Dev.expl[i] <- summary(fit)$dev.expl
 z <- range(fitted(fit))
 if(z[2] > abs(1.75 * z[1]))  ZLIM <- c(-4/7*z[2], z[2]) else ZLIM <- c(z[1], -7/4 * z[1]) 
 vis.gam(fit, c("sp.long", "sp.lat"), n.grid=100, plot.type="contour", color="topo", nlevels=8, labcex=0.8,
   too.far = 0.04, xlab = "Longitude", ylab = "Latitude", xlim=c(-178,-157), ylim=c(54,61),
   zlim = ZLIM, main = "", cex.lab=1.4, cex.main=1.6)
 map('worldHires',fill=T,xlim=c(-178,-157), ylim=c(54,61),add=T)
 text(-159, 60.2, "Alaska", cex=1.5, col = "white")
 lines(EBS.contours$z50, lty=3, lwd=2)
 lines(EBS.contours$z100, lty=2, lwd=2)
 lines(EBS.contours$z200, lty=1, lwd=2)

# Plot time series of first PC
 x <- barplot(pc1)
 axis(1, at=x, labels=years, las=2, col=0, line=-0.5)
 PC1[,i] <- pc1
 par(usr = c(par("usr")[1:2], range(cpue.ann[,i])))
 lines(x, cpue.ann[,i], col=4)        # Superimpose average CPUE by year!
 par(usr = c(par("usr")[1:2], range(avg.T[,i])))
 lines(x, avg.T[,i], col=2)           # Superimpose average temperature by year!

# Species:
 mtext(paste(Taxa[i,2], " (", Taxa[i,1], ")", sep=""), 3, outer=T, cex=2)
}

dev.off()     #!!!!!! make sure to run to redirect output back to graphics window

# clean up a bit:
rm(anom, cpue, x, K, i, j, ev, f, cf, n, sp.lat, sp.long)

# Output from 'cgi'
names(CGI)



############################################################################
############################################################################

### Not sure if stuff below still works but not essential!

# I did PCA of total CPUE and of BT as well as PCA by stratum

############################################################################
# Total CPUE:

win.graph(9,7.5)
# To save graphical output to pdf file:
# pdf("TotalCPUE PCA.pdf", width=9, height=7.5)

par(omi=c(0,0,0.5,0), mar=c(4.2,4.2,2,1))
layout(matrix(c(2,3,0,1),2,2), width = c(2,1), height=c(1.5,1))

cpue <- BShauls$Total.CPUE
cpue <- tapply(cpue, list(BShauls$YEAR, BShauls$STATION), mean, na.rm=T)
cpue <- cpue[, st]  # Select stations sampled each year
if(sum(is.na(cpue)) > 0) stop("Missing station/year combination in cpue")

# Compute and save diversity of catches across stations: 
 f <- function(x) {
  x <- x[x > 0]
  pi <- x/sum(x)
  - sum(pi * logb(pi, base = 2))
 }
 Div.CPUE <- apply(cpue,1,f) 

 j <- is.element(BShauls$STATION, dimnames(cpue)[[2]])
 Avg.T <- tapply(BShauls$GEAR_TEMP[j], BShauls$YEAR[j], mean, na.rm=T)

# Remove annual means across stations (CPUE anomalies)
 anom <- sweep(cpue^0.25, 1, apply(cpue^0.25,1,mean))
 anom.svd <- svd(var(anom))
 ev <- anom.svd$d
 evec <- anom.svd$u[,1]
 x <- barplot(ev[1:8]/sum(ev)*100, xlab = "Principal Component", ylab="Proportion of variance")
 axis(1, at=x, 1:8, col=0)
 title("Scree plot")
 ev[1:5]/sum(ev) * 100  # Variance explained (%)

# Spatial pattern:
 fit <- gam(evec ~ s(lat, long, k=60))
 summary(fit)
 z <- range(fitted(fit))
 if(z[2] > abs(1.75 * z[1]))  ZLIM <- c(-4/7*z[2], z[2]) else ZLIM <- c(z[1], -7/4 * z[1]) 
 vis.gam(fit, c("long", "lat"), n.grid=100, plot.type="contour", color="topo", nlevels=8, labcex=0.8,
   too.far = 0.04, xlab = "Longitude", ylab = "Latitude", xlim=c(-178,-157), ylim=c(54,61),
   zlim = ZLIM, main = "", cex.lab=1.4, cex.main=1.6)
 map('worldHires',fill=T,xlim=c(-178,-157), ylim=c(54,61),add=T)
 text(-159, 60.2, "Alaska", cex=1.5, col = "white")
 lines(EBS.contours$z50, lty=3, lwd=2)
 lines(EBS.contours$z100, lty=2, lwd=2)
 lines(EBS.contours$z200, lty=1, lwd=2)

# Time series of first PC:
 pc1 <- as.vector(scale(anom %*% evec))
 x <- barplot(pc1)
 axis(1, at=x, labels=years, las=2, col=0, line=-0.5)
 mtext("PC 1", side=2, cex=1.6, line=2.5)
 pc1
 par(usr = c(par("usr")[1:2], range(Avg.CPUE)))
 lines(x, Avg.CPUE, lty=2)
 par(usr = c(par("usr")[1:2], range(Avg.T)))
 lines(x, Avg.T,lty=2, col=2) 

cor(cbind(Avg.CPUE, Avg.T, pc1))
cor.test(Avg.T, pc1)

 mtext("Total CPUE", outer=T, cex=2)


############################################################################
# PCA of temperature:

bt <- BShauls$GEAR_TEMP
bt <- tapply(bt, list(BShauls$YEAR, BShauls$STATION), mean, na.rm=T)
bt <- bt[, st]  # Select stations sampled each year

# Compute and save average bottom temperature (BT)
 Avg.BT <- apply(bt,1,mean, na.rm=T)

# Remove annual means across stations (to get anomalies relative to annual mean)
 anom <- sweep(bt, 1, apply(bt,1,mean, na.rm=T))
 anom.svd <- svd(cov(anom, use="pair"))
 ev <- anom.svd$d
 evec <- anom.svd$u[,1]
 x <- barplot(ev[1:8]/sum(ev)*100, xlab = "Principal Component", ylab="Proportion of variance")
 axis(1, at=x, 1:8, col=0)
 title("Scree plot")
 ev[1:5]/sum(ev) * 100  # Variance explained (%)

# Spatial pattern:
 fit <- gam(evec ~ s(lat, long, k=60))
 summary(fit)
 z <- range(fitted(fit))
 if(z[2] > abs(1.75 * z[1]))  ZLIM <- c(-4/7*z[2], z[2]) else ZLIM <- c(z[1], -7/4 * z[1]) 
 vis.gam(fit, c("long", "lat"), n.grid=100, plot.type="contour", color="topo", nlevels=8, labcex=0.8,
   too.far = 0.04, xlab = "Longitude", ylab = "Latitude", xlim=c(-178,-157), ylim=c(54,61),
   zlim = ZLIM, main = "", cex.lab=1.4, cex.main=1.6)
 map('worldHires',fill=T,xlim=c(-178,-157), ylim=c(54,61),add=T)
 text(-159, 60.2, "Alaska", cex=1.5, col = "white")
 lines(EBS.contours$z50, lty=3, lwd=2)
 lines(EBS.contours$z100, lty=2, lwd=2)
 lines(EBS.contours$z200, lty=1, lwd=2)

# Time series of first PC (ignore missing values):
 pc1 <- t(anom)
 pc1 <- pc1 * evec
 pc1 <- as.vector(scale(apply(pc1, 2, sum, na.rm=T)))
 x <- barplot(pc1)
 axis(1, at=x, labels=years, las=2, col=0, line=-0.5)
 pc1
 par(usr = c(par("usr")[1:2], range(Avg.BT)))
 lines(x, Avg.BT, lty=2)
 
 mtext("Bottom temperature anomalies", outer=T, cex=2)

 BT.pc1 <- pc1

############################################################################
# Analysis by stratum:

table(BShauls$STRATUM, BShauls$YEAR)
LAT <- tapply(BShauls$LATITUDE, BShauls$STRATUM, mean, na.rm=T)[1:10]
LONG <- tapply(BShauls$LONGITUDE, BShauls$STRATUM, mean, na.rm=T)[1:10]
map.EBS()
text(LONG, LAT, names(LAT), cex=2)
map.EBS(img=F)
x <- BShauls[BShauls$YEAR == 2009 & BShauls$STRATUM < 70, c("LONGITUDE", "LATITUDE", "STRATUM")]
points(x$LONG, x$LAT, pch=16, col=factor(x$STRATUM))

cpue <- BShauls$Total.CPUE
cpue <- tapply(cpue, list(BShauls$YEAR, BShauls$STRATUM), mean, na.rm=T)
cpue <- cpue[, as.numeric(dimnames(cpue)[[2]]) < 70]  # Select strata sampled each year

cpue.pca <- prcomp(cpue, center=T, scale=T)
summary(cpue.pca)
plot(cpue.pca)
cpue.pca$rotat[,1:3]

# Correlation between average bottom temperature and 1st PC:
plot(years, -scale(cpue.pca$x[,1]), type="b")
lines(years, scale(Avg.BT),col=2)
cor(Avg.BT, cpue.pca$x[,1])^2


