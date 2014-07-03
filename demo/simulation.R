params000 <- list(n=c(100,200,300),
                  beta0=5,
                  beta1=10,
                  beta2=1,
                  sigmaV=0.5, 
                  sigmaU=2.5)

params000T <- c(params000, list(mu=1))
params100 <- c(params000, list(rhoY=0.2))
params100T <- c(params000T, list(rhoY=0.2))

params110 <- c(params100, list(rhoV=0.2))
params101 <- c(params100, list(rhoU=0.2))
params111 <- c(params110, list(rhoU=0.2))

params010 <- params110
params010$rhoY <- NULL

params011 <- params111
params011$rhoY <- NULL

params001 <- params011
params001$rhoV <- NULL

ctrl <- list(true.initial=F, seed=999, cores=detectCores())
res000 <- ezsimspfrontier(100, params = params000,  inefficiency = "half-normal",logging = "info", control=ctrl)
save(res000, file="res000.rData")

res000T <- ezsimspfrontier(100, params = params000T, inefficiency = "truncated",logging = "info", control=ctrl)
save(res000T, file="res000T.rData")

res100 <- ezsimspfrontier(100, params = params100,  inefficiency = "half-normal",logging = "info", control=ctrl)
save(res100, file="res100.rData")

res100A <- ezsimspfrontier(100, params = params100,  inefficiency = "half-normal",logging = "info", 
                               control=c(ctrl,list(ignoreWy=T)))
save(res100A, file="res100A.rData")

res100T <- ezsimspfrontier(100, params = params100T, inefficiency = "truncated",logging = "info", control=ctrl)
save(res100T, file="res100T.rData")

params001A <- c(params001, list(rhoY=0))
res001A <- ezsimspfrontier(100, params = params001A, inefficiency = "half-normal",logging = "info", 
                           control=c(ctrl,list(ignoreWu=T,replaceWyWu=T)))
save(res001A, file="res001A.rData")

params010A <- c(params010, list(rhoY=0))
res010A <- ezsimspfrontier(100, params = params010A, inefficiency = "half-normal",logging = "info", 
                           control=c(ctrl,list(ignoreWv=T,replaceWyWv=T)))
save(res010A, file="res010A.rData")



#All tests above work as appropriate

res010 <- ezsimspfrontier(100, params = params010, inefficiency = "half-normal",logging = "info",control=ctrl)
res001 <- ezsimspfrontier(100, params = params001, inefficiency = "half-normal",logging = "info",control=ctrl)
res101 <- ezsimspfrontier(100, params = params101, inefficiency = "half-normal",logging = "info",control=ctrl)
res111 <- ezsimspfrontier(100, params = params111, inefficiency = "half-normal",logging = "info",control=ctrl)



ctrl <- list(true.initial=TRUE, seed=0, cores=detectCores())
params001$n <- c(100)
res001.100 <- ezsimspfrontier(100, params = params001, inefficiency = "half-normal",logging = "info",control=ctrl)
save(res001.100, file="res001.100.rData")

params001$n <- c(200)
res001.200 <- ezsimspfrontier(100, params = params001, inefficiency = "half-normal",logging = "info",control=ctrl)
save(res001.200, file="res001.200.rData")

params001$n <- c(300)
res001.300 <- ezsimspfrontier(100, params = params001, inefficiency = "half-normal",logging = "info",control=ctrl)
save(res001.300, file="res001.300.rData")


