require(spfrontier)

data( airports )
W <- constructW(cbind(airports$lon, airports$lat),airports$ICAO_code)

formula <- log(PAX) ~ log(runways) + log(checkins) +log (gates)
ols <- lm(formula , data=airports)
summary(ols )
plot(density(stats::residuals(ols)))
skewness(stats::residuals(ols))
 
model <- spfrontier(formula , data=airports)
summary(model )

model <- spfrontier(formula , data=airports, W_y=W)
summary(model )

# Takes a long time
#model <- spfrontier(formula , data=airports, W_y=W, W_v=W, logging="debug",control=list())
#summary(model )
