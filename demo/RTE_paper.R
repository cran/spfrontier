rm(list=ls())

library(spdep)
library(pastecs)
library(car)
library(maptools)
library(maps)
library(ggmap)
library(corpcor)
library(spgwr)
library(raster)
library(RColorBrewer)
library(classInt)
library(frontier)
library(spfrontier)
library(gstat)



rowStdrt = function(W){
    for (j in 1:nrow(W)){
        W[j,] <- W[j,]/sum(W[j,])
    }
    return(W)
}


constructW <- function(coords, labels){
    W <- 1/as.matrix(dist(coords))
    colnames(W) <- labels
    rownames(W) <- labels
    W[which(W==Inf)] <- 0
    return(W)
}
#data <-read.csv("data/data.csv", header=T)
#save(data, file = "data/RTEpaper.rda")
data(RTEpaper, envir = environment())
data<-RTEpaper

#data <- subset(data, PAX > 4000000)





data$PAXperRunway <- data$PAX/(data$RunwayCount*1000000)
data$PAXperRunwayLength <- data$PAX/(data$RunwayLength*1000)
data$PAXperPop100 <- data$PAX/(data$Population100km)
data$TouristsPerkm<- data$Tourists/data$Area
#data <- subset(data, data$ICAO!="EGLC")
#data <- subset(data, data$ICAO!="ELLX")
stat.desc(data)

nrow(data)

data$OwnPublic <- ifelse(data$ownership=="public",1,0)
data$OwnMinorPrivate <- ifelse(data$ownership=="minor private",1,0)
data$OwnMajorPrivate <- ifelse(data$ownership=="major private",1,0)
data$OwnPrivate <- ifelse(data$ownership=="private",1,0)



W <- constructW(cbind(data$longitude, data$latitude),data$ICAO)

W <- rowStdrt(W)
listw <- mat2listw(W)

#formula <- log(PAX) ~ log(RunwayCount) +  log(Population100km)+  log(GDPpc)+  log(TouristsPerkm) + SouthIsland + hub + international
#formula <- log(PAX) ~ log(RunwayLength) +  log(Population100km)+  log(GDPpc)+  log(TouristsPerkm) + hub + international

#formula <- log(PAX) ~ log(Population100km) + log(Routes)+ log(RunwayCount) + log(CheckinCount)+ log(GateCount)+ log(GDPpc) + Island 

formula <- -log(PAX) ~ log(0.00001+cargo/PAX) + log(RunwayLength) +  log(Population100km)+  log(GDPpc)+  log(TouristsPerkm) + hub + international

eu.sfa <- spfrontier( formula , data=data, logging="info", costFrontier=T)
summary(eu.sfa)

ini <- c(coefficients(eu.sfa)$beta,0,coefficients(eu.sfa)$sigmaV, coefficients(eu.sfa)$sigmaU,0)

eu.ssf1000 <- spfrontier( formula , W_y=W,data=data, logging="info", costFrontier=T)
summary(eu.ssf1000)

eu.ols <- lm( formula , data=data)
summary(eu.ols)
vif(eu.ols)

resid <- stats::residuals(eu.ols)
dens <- density(resid)
dens <-  density(spfrontier::efficiencies(eu.sfa))
dens <-  density(spfrontier::efficiencies(eu.ssfa0010 ))


plot(dens,main="OLS residuals' kernel density")
x1 <- min(dens$x)  
x2 <- max(dens$x)

dd <- with(dens,data.frame(x,y))
library(ggplot2)
qplot(x,y,data=dd,geom="line",ylab="density",xlab="OLS residuals")+
  geom_ribbon(data=subset(dd,x>x1 & x<x2),aes(ymax=y),ymin=0,
              fill="red",colour=NA,alpha=0.2)

hist(spfrontier::efficiencies(eu.ssfa0010 ))
skewness(resid)



eu.sfa <- sfa( formula , data=data, ineffDecrease=FALSE)
summary(eu.sfa)


moran.test(resid,listw,randomisation=FALSE,alternative="two.sided")
moran.plot(resid,listw)
moran.plot(as.vector(spfrontier::efficiencies(eu.sfa)),listw)

moran.test(spfrontier::efficiencies(eu.sfa),listw,randomisation=FALSE,alternative="two.sided")

lm.LMtests(eu.ols, listw, test=c("LMerr","RLMerr","LMlag","RLMlag","SARMA"))

summary(lagsarlm(formula ,data=data, listw))
eu.sem <- errorsarlm(formula ,data=data, listw)
summary(eu.sem)

resid <- residuals(eu.sem) - eu.sem$lambda * W %*% residuals(eu.sem)

eu.sfa <- spfrontier( formula , data=data, logging="info", costFrontier=T)
ini <- c(coefficients(eu.sfa)$beta,coefficients(eu.sfa)$sigmaV, coefficients(eu.sfa)$sigmaU)
eu.ssfa0010 <- spfrontier( formula , data=data,W_v=W, logging="debug", 
                           initialValues=c(ini,0.5)+0.01, costFrontier=T,
                           control=list(optim.control = list(tol=1e-2, iterlim=2000,repeatNM=1,printLevel=999)))
summary(eu.ssfa0010)
save(eu.ssfa0010, file="eu.ssfa0010-review2-2.rData")
load(file="eu.ssfa0010-review.rData")

ini <- c(coefficients(eu.sfa)$beta, coefficients(eu.sfa)$sigmaV, coefficients(eu.sfa)$sigmaU)
eu.ssfa0001 <- spfrontier( formula , data=data, W_u=W, logging="debug", 
                           initialValues=c(ini, 0),
                           control=list(optim.control = list(tol=1e-2, iterlim=2000,repeatNM=1)))
summary(eu.ssfa0001)
save(eu.ssfa0001, file="eu.ssfa0001.rData")


eu.ssfa1000 <- spfrontier( formula , data=data, W_y=W)
summary(eu.ssfa1000)

eu.bw <- gwr.sel(formula, data=data,coords=cbind(data$latitude, data$longitude))

eu.gwr <- gwr(formula, data=data,coords=cbind(data$latitude, data$longitude), bandwidth=eu.bw, hatmatrix=TRUE)
eu.gwr
eu.gwr$SDF$latitude <- data$latitude
eu.gwr$SDF$longitude <- data$longitude
eu.gwr$SDF$ICAO <- data$ICAO
eu.gwr$SDF$log_RunwayLength_t <- eu.gwr$SDF[["log(RunwayLength)"]]/eu.gwr$SDF[["log(RunwayLength)_se"]]
eu.gwr$SDF$log_Population100km_t <- eu.gwr$SDF[["log(Population100km)"]]/eu.gwr$SDF[["log(Population100km)_se"]]
eu.gwr$SDF$log_GDPpc_t <- eu.gwr$SDF[["log(GDPpc)"]]/eu.gwr$SDF[["log(GDPpc)_se"]]
eu.gwr$SDF$log_TouristsPerkm_t <- eu.gwr$SDF[["log(TouristsPerkm)"]]/eu.gwr$SDF[["log(TouristsPerkm)_se"]]

spplot(eu.gwr$SDF, "gwr.e", add=TRUE)
spplot(eu.gwr$SDF, "log.RunwayCount.")
spplot(eu.gwr$SDF, "log_RunwayCount_t")
spplot(eu.gwr$SDF, "log_Population100km_t")
spplot(eu.gwr$SDF, "log_GDPpc_t")
spplot(eu.gwr$SDF, "log.Population100km.")
spplot(eu.gwr$SDF, "log.GDPpc._se")


wd <- "D:/Dmitry/Dropbox/science/PhD/R/spfrontier/phd"
setwd(wd)
source("utils.R")

library(maptools)
library(classI)

world<-readShapeSpatial("TM_WORLD_BORDERS-0.3/TM_WORLD_BORDERS-0.3.shp", proj4string=CRS("+proj=longlat"))


world <- world[!is.na(world$FIPS),]

europe <- world[world$REGION == 150 | world$FIPS == "CY",]
eu <- europe

nclr <- 10
plotclr <- colorRampPalette(c("yellow", "red"))(nclr)
classI <- classIntervals(eu.gwr$SDF[["log_RunwayCount_t"]], nclr, style = "equal")
brks <- round(classI$brks, 1)
classI <- classIntervals(eu.gwr$SDF[["log_RunwayCount_t"]], nclr, style = "fixed",fixedBreaks = brks)
colcode <- findColours(classI, plotclr)

plot(eu, border="grey", xlim = c(0 , 10) , ylim = c(30 , 70))
points(eu.gwr$SDF$longitude,eu.gwr$SDF$latitude,pch=20,col=colcode)
legend("bottomleft", legend = names(attr(colcode, "table")), title="Title",fill = attr(colcode, "palette"), cex = 1, bty = "n")


nclr <- 7
dat <- residuals(eu.ols)
dat <- spfrontier::efficiencies(eu.sfa)
dat <- eu.gwr$SDF[["log(TouristsPerkm)"]]
dat <- eu.gwr$SDF[["log(Population100km)"]]
dat <- eu.gwr$SDF[["log(GDPpc)"]]
dat <- eu.gwr$SDF[["log(RunwayLength)"]]

dat <- eu.gwr$SDF$log_TouristsPerkm_t
dat <- eu.gwr$SDF$log_Population100km_t
dat <- eu.gwr$SDF$log_GDPpc_t
dat <- eu.gwr$SDF$log_RunwayLength_t
dat <- data$effSSFA-data$effSFA

par( mfrow = c(1, 1 ) )

plotclr <- colorRampPalette(c("yellow", "red"))(nclr ) 
classI <- classIntervals(dat, nclr, style = "equal")
brks <- round(classI$brks, 1)
classI <- classIntervals(dat, nclr, style = "fixed",fixedBreaks = brks)
colcode <- findColours(classI, plotclr)
plot(eu, border="darkgrey",col="gray90", xlim = c(-5 , 0) , ylim = c(29 , 70))
points(data$longitude,data$latitude,pch=21,col="darkgrey",bg=colcode,cex=1)
legend("bottomleft", legend = names(attr(colcode, "table")), title="Tourists: t-values",fill = attr(colcode, "palette"), cex = 1, bty = "n")



data$effSFA <- efficiencies(eu.sfa)
data$effSSFA <- spfrontier::efficiencies(eu.ssfa0010)
boxplot(as.vector(data$effSFA), as.vector(data$effSSFA),names=c("SFA","SSFA"),main="Efficiency levels")

d<-data[order(data$effSSFA-data$effSFA , decreasing=T),]
head(subset(d, select=c(ICAO,AirportName,Country,PAX,effSFA,effSSFA)),5)
tail(subset(d, select=c(Country,ICAO,AirportName,PAX,effSFA,effSSFA)),5)

d<-data[order(data$effSSFA , decreasing=T),]
subset(d, select=c(ICAO,effSFA,effSSFA))

d<-data[order(data$PAXperRunwayLength , decreasing=T),]
subset(d, select=c(ICAO,PAXperRunwayLength))


coordinates(data) =  ~longitude + latitude
proj4string(data) =  "+proj=longlat +datum=WGS84" 

x.range <- as.numeric(c(min(data$longitude), max(data$longitude))) 
y.range <- as.numeric(c(min(data$latitude), max(data$latitude))) 
grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = 0.1), y = seq(from = y.range[1], to = y.range[2], by = 0.1))  
coordinates(grd) <- ~ x + y
proj4string(grd) =  "+proj=longlat +datum=WGS84" 
gridded(grd) <- TRUE
plot(grd, cex = 1.5, col = "grey")
points(data, pch = 1, col = "red", cex = 1)
idw1 <- idw(formula = data$effSFA ~ 1, locations = data, newdata = grd,idp=2)
idw.output <- as.data.frame(idw1)
names(idw.output)[1:3]<- c("long", "lat", "var1.pred")
dat <- data$effSFA
nclr <- 7
plotclr <- colorRampPalette(c("cyan", "orange"))(nclr ) 
classI <- classIntervals(dat, nclr, style = "equal")
brks <- round(classI$brks, 1)
classI <- classIntervals(dat, nclr, style = "fixed",fixedBreaks = brks)
colcode <- findColours(classI, plotclr)
ggplot() + geom_tile(data = idw.output, alpha = 1, aes(x = long, y = lat, 
    fill = round(var1.pred, 1))) + scale_fill_gradient(low = "cyan", high = "orange") + 
    geom_path(data = eu, aes(long, lat, group = group), colour = "black") + xlim(x.range) +ylim(y.range)+
    geom_point(data = as.data.frame(data), aes(x = longitude, y = latitude), shape = 21, fill= colcode, col= colcode)+
	labs(fill = "Efficiency level", title = "Efficiency of EU airports (SF model)")+ theme_bw()




v1 <- variogram(data$effSFA~1, data, cutoff=500, width=10)
lo <- loess(v1$gamma~v1$dist, span=1.5)
plot(v1$dist,v1$gamma, xlab="Distance", ylab="Semi-variance",pch=20, col="blue")
lines(v1$dist,predict(lo), col='red', lwd=2)

ggplot() + geom_point(data = as.data.frame(v1), aes(x = dist, y = gamma), shape = 20, fill= "blue", col="blue")+
geom_line(data = as.data.frame(v1), aes(x = dist, y = predict(lo)), col="red", size=2) +
labs(x = "Distance",y="Semi-variance", title = "Variogram of SF efficiencies values")


geo.dists <- dist(cbind(data$longitude, data$latitude))
value.dists <- dist(data$effSFA)
mantel <- mantel.rtest(geo.dists, value.dists, nrepet = 999)




