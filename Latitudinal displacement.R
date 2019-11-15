# Average center of gravity (latitude) across 46 taxa

CG.lat
x <- scale(CG.lat, scale=T)
y <- 1982:2012
plot(y, x <- apply(x,1,mean,na.rm=T), type="b")

par(mar=c(3,5,1,1))
plot(y,x, type="b", col=4, xlab = "", ylab = "Normalized northward shift", 
   cex=1.5, pch=16, cex.lab = 1.8, cex.axis=1.4)
abline(fit<-gls(x ~ y, correlation=corAR1()), col=2, lwd=3)
summary(fit)

map.EBS()
j <- BShauls$YEAR == 2008
points(BShauls$LONG[j], BShauls$LAT[j], pch=16, cex=0.7) #, col="yellow")


