# Copyright 2014 by Dmitry Pavlyuk <Dmitry.V.Pavlyuk@gmail.com>

#
# Routines for maximum likelihood estimator
#

#
# Wrapper for 'optim' procedure from 'stats' package
#
optimEstimator <- function(formula, data, func.lf, ini,gr=NULL){
    ret <- prepareEstimates(status = 1000)
    if (is.null(ini)){
        logging("Ini is not defined",caller = match.call())
        return(prepareEstimates(status = 1001))
    }
    tryCatch({
        control <- envirGet("control")
        #constr <- envirGet("constr")
        pscale <- envirGet("parscale")
        lowerBounds <- envirGet("lowerBounds")
        upperBounds <- envirGet("upperBounds")
        if (is.null(lowerBounds)){
            lowerBounds<- -Inf
        }
        if (is.null(upperBounds)){
            upperBounds<- Inf
        }
        control$optim.control$maximize <- TRUE
        control$optim.control$kkt <- (control$hessian=="optim")
        control$optim.control$starttests <- FALSE
        #control$optim.control$parscale <- pscale
        if (is.null(gr)){
            method <- "Nelder-Mead"
            lowerBounds <- -Inf
            upperBounds <- Inf
        }else{
            method <- "L-BFGS-B"
            control$optim.control$reltol <- NULL
            control$optim.control$abstol <- NULL
        }
            p<-optimx(par=ini,fn=func.lf, gr=gr,
                      lower = lowerBounds,
                      upper = upperBounds,
                      hessian = (control$hessian=="optim"),
                      method = method,
                    #constraints=constr,
                 #reltol=control$optim.control$tol,
                 #printLevel=control$optim.control$printLevel,
                 #parscale=pscale,
                 #iterlim=control$optim.control$iterlim
                 control=control$optim.control
                 )
            if (p$convcode==0){
                hessian <- NULL
                if (control$hessian=="optim"){
                    hessian <- attr(p, "details")[ ,"nhatend"][[1]]
                }else if (control$hessian=="numeric"){
                        N <- length(coef(p))
                        calls <- 2 + 4 * (N^2 + N)
                        logging(paste("Calculating numeric hessian, number of required function calls is ",calls," ..."), level = "info")
                        hessian <- numDeriv::hessian(func.lf, coef(p))
                    }
                attr(p, "calchessian") <- hessian
            }
        #p <- psoptim(par=NA, fn=func.lf, lower=c(-Inf, -Inf, -Inf, 0,0,-1), upper=c(Inf, Inf, Inf, Inf,Inf,1))
        logging(summary(p), level = "debug")
        logging(paste("Convergence is",
                      ifelse(p$convcode==0, "", "NOT"),"achieved",
                      "[",p$convcode,"]"),
                level = "info",
                caller = match.call())    
        logging(paste("Log-Likelihood value = ",p$value),
                level = "info",
                caller = match.call())  
        ret <- prepareEstimates(estimates = p)
    }, error = function(e){ 
        logging(e$message, level="warn")
    })
    return(ret)
}

#succcessCode <- function(p){
#    res <- p$code
#    if (p$type=="Newton-Raphson maximisation"){
#        if (res<3) res<-0
#    }
#    return(res)
#}
#
# Creates the ModelEstimates class on the base of raw (optim) estimation results
#
prepareEstimates <- function(estimates = NULL, status = 0){
    ret <- new("ModelEstimates", status = status)
    if (!is.null(estimates)) {
        if (estimates$value == Infin){
            return(new("ModelEstimates", status = 1002))
        }
        ret <- new("ModelEstimates",
                   resultParams = coef(estimates),
                   status = estimates$convcode,
                   logL = estimates$value,
                   logLcalls = c(estimates$fevals),
                   hessian = attr(estimates, "calchessian") 
        )
        #Setting hessian and implicitly calculating standard errors
        #if (!is.null(estimates$hessian)){
        #    hessian(ret) <- estimates$hessian
        #}
    }
    return(ret)
}