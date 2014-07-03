library(spdep)
library(car)

rm(list=ls())
data(airports, envir = environment())

########################
#ALL 

airports2011 <- subset(airports, Year==2011)
data <- airports2011 
nrow(data)
W <- constructW(cbind(data$longitude, data$latitude),data$ICAO)

formula <- log(PAX) ~ log(Population100km) + log(Routes) + log(GDPpc) +Island
formula <- log(PAX+0.0001) ~ log(Cargo/PAX+0.0001) + log(Population100km) + log(Routes) + log(GDPpc)


##########################

########################
#SPAIN 

airports2010.Spain <- subset(airports, Year==2010 & Country=="Spain")
nrow(airports2010.Spain)
airports2010.Spain.complete <- with(airports2010.Spain, airports2010.Spain[complete.cases(RevenueTotal,ATM,PAX,DA,StaffCost,TerminalCount,RunwayCount,Population100km,Routes,longitude,latitude),])
data <- airports2010.Spain.complete 
nrow(data)
W <- constructW(cbind(data$longitude, data$latitude),data$ICAO)

formula <- log(RevenueTotal) ~ log(ATM) + log(PAX) + log(DA) + log(StaffCost)+ log(TerminalCount)+ log(RunwayCount) + log(Population100km) + log(Routes)
formula <- log(PAX+0.0001) ~ log(Routes) + log(TerminalCount) + log(RunwayCount)+ log(Population100km) + Island

formula <- log(RevenueTotal) ~ log(ATM) + log(PAX) + log(StaffCost)+ log(TerminalCount)+ log(Population100km) + log(Routes)
#Good!


##########################

########################
#UK 

airports2012.UK <- subset(airports, Year==2012 & Country=="United Kingdom")
nrow(airports2012.UK)
data <- airports2012.UK
nrow(data)

W <- constructW(cbind(data$longitude, data$latitude),data$ICAO)

formula <- log(PAX+0.0001) ~ log(Routes) + log(Population100km) + Island


##########################

########################
#Greece 

rm(list=ls())
data(airports.greece, envir = environment())
data <- airports.greece


W <- constructW(cbind(data$lon, data$lat),data$ICAO)

formula <- log(WLU) ~ log(openning_hours) + log(runway_area) + log(terminal_area) + log(parking_area) + island + international
#Good!

formula <- log(WLU) ~ log(openning_hours) + log(terminal_area)

##########################

ols <- lm(formula , data=data)
summary(ols )
vif(ols)
resid <- stats::residuals(ols)
plot(density(resid))
skewness(resid)

listw <- mat2listw(W)

moran.test(resid,listw,randomisation=FALSE,alternative="two.sided")
lm.LMtests(ols, listw, test=c("LMerr","RLMerr","LMlag","RLMlag","SARMA"))

summary(lagsarlm(formula ,data=data, listw))
summary(errorsarlm(formula ,data=data, listw))

 
model000 <- spfrontier( formula , data=data)
summary(model000 )

model100 <- spfrontier(formula , data=data, W_y=W)
summary(model100 )

model010 <- spfrontier(formula , data=data, W_v=W, logging="debug",control=list())
summary(model010)

model001 <- spfrontier(formula , data=data, W_u=W, logging="debug",control=list())
summary(model001)