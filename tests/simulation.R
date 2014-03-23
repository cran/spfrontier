require(spfrontier)

beta0 = beta1 = beta2 = sigmaV = sigmaU = n = sigmaX = rhoY = rhoV = rhoU = mu = NULL

params000 <- list(n=c(100),
                  sigmaX=10, 
                  beta0=1,
                  beta1=-2,
                  beta2=3, 
                  sigmaV=0.2, 
                  sigmaU=0.75)

params000T <- params000
params000T$Mu <- 0.4

params100 <- params000
params100$rhoY <- 0.6

params100T <- params000T
params100T$rhoY <- 0.6


params110 <- params100
params110$rhoV <- 0.7

params010 <- params110
params010$rhoY <- NULL

params111 <- params110
params111$rhoU <- 0.5

params011 <- params111
params011$rhoY <- NULL

params001 <- params011
params001$rhoV <- NULL

#res <- ezsimspfrontier(100, params = params000,  seed = 999, inefficiency = "half-normal",logging = "info")
#res <- ezsimspfrontier(100, params = params000T, seed = 999, inefficiency = "truncated",logging = "info")
#res <- ezsimspfrontier(100, params = params100,  seed = 999, inefficiency = "half-normal",logging = "info")
#res <- ezsimspfrontier(100, params = params100T, seed = 999, inefficiency = "truncated",logging = "info")
#All tests above work as appropriate
#res <- ezsimspfrontier(10, params = params010, seed = 999, inefficiency = "half-normal",logging = "debug")
#A problem with sigmaV
#res <- ezsimspfrontier(10, params = params001, seed = 999, inefficiency = "half-normal",logging = "info")
#res <- ezsimspfrontier(10, params = params110, seed = 999, inefficiency = "half-normal",logging = "info")
#res <- ezsimspfrontier(10, params = params111, seed = 999, inefficiency = "half-normal",logging = "info")
#res <- ezsimspfrontier(10, params = params011, seed = 999, inefficiency = "half-normal",logging = "info")