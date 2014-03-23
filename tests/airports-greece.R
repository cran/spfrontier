require(spfrontier)

data(airports.greece)
formula <- log(WLU) ~ log(openning_hours) + log(runway_area) + log(terminal_area)
W <- constructW(cbind(airports.greece$lon, airports.greece$lat),airports.greece$ICAO)

model000 <- spfrontier(formula , data=airports.greece, logging="info")
summary(model000 )
eff000 <- efficiencies(model000)

model100 <- spfrontier(formula , data=airports.greece, logging="info", W_y=W)
summary(model100 )
eff100 <- efficiencies(model100)

# Takes a long time
#model010 <- spfrontier(formula , data=airports.greece, logging="debug", W_v=W)
#summary(model010 )
#eff010 <- efficiencies(model010)

# Takes a long time
#model110 <- spfrontier(formula , data=airports.greece, logging="debug", W_y=W, W_v=W)
#summary(model110)
#eff110 <- efficiencies(model110)

# Takes a long time
#model001 <- spfrontier(formula , data=airports.greece, logging="debug", W_u=W)
#summary(model001)
#eff001 <-efficiencies(model001)

#cbind(eff000, eff100, eff010, eff110, eff001)