beta0 = beta1 = beta2 = sigmaV = sigmaU = n = sigmaX = rhoY = rhoV = rhoU = mu =loggingLevel= inefficiency=parDef=NULL

spfrontier.dgp <- function(){
    if (!is.null(mu)){
        if (abs(mu)>sigmaU){
            cat("DGP: Truncated normal mean (mu) is higher than standard deviation, which can lead to non-skewed residuals")
        }
    }
    
    formula <- as.formula("y ~ X1 + X2")
    beta<-c(beta1,beta2)
    k <- length(beta)
    X <- matrix(rnorm(n*k,0,sigmaX),n, k)
    
    W_v <- NULL
    SpW2 <- diag(n)
    if (!is.null(rhoV)){
        W_v <- genW(n,type="rook")
        SpW2 <- solve(diag(n)-rhoV*W_v)
    }
    v <-  SpW2%*%rmvnorm(1,mean = rep(0, n),sigma = sigmaV^2*diag(n))[1,]
    if (!is.null(rhoV)){
        print(coef(lm(v~Wv-1, data=data.frame(v, Wv = W_v%*%v))))
    }
    
    W_u <- NULL
    muval <- 0
    if (!is.null(mu))
        muval <- mu
    SpW3 <- diag(n)
    if (!is.null(rhoU)){
        W_u <- genW(n,type="queen")
        SpW3 <- solve(diag(n)-rhoU*W_u)
    }
    u <- SpW3%*%rtmvnorm(1,mean = rep(muval, n),sigma = sigmaU^2*diag(n),algorithm="gibbs", lower=rep(0, n))[1,]
    if (!is.null(rhoU)){
        print(coef(lm(u~Wu-1, data=data.frame(u, Wu = W_u%*%u))))
    }
    sk <- skewness(v-u)
    if (sk>=0){
        cat("DGP: Skewness of generated residuals is non-negative = ",sk)
    }
    
    #plot(density(v - u))
    y <- beta0 + X %*% beta + v - u
    
    W_y <- NULL
    if (!is.null(rhoY)){
        W_y <- genW(n,type="queen")
        SpW <- solve(diag(n)-rhoY*W_y)
        y <- SpW%*%y
    }
    dat <- data.frame(y,X)
    colnames(dat) <-c('y',paste("X", seq(k), sep = ""))
    tv <- evalFunctionOnParameterDef(parDef,spfrontier.true.value)
    result <- list(formula=formula, data=dat,W_y=W_y,W_v=W_v,W_u=W_u, tv=tv,
                   loggingLevel=loggingLevel,inefficiency=inefficiency)
    return(result)
}

spfrontier.estimator <- function(d){
    modelEstimates <- spfrontier(d$formula,d$data,W_y=d$W_y,W_v=d$W_v,W_u=d$W_u,
                                 logging = d$loggingLevel,inefficiency=d$inefficiency,onlyCoef=T,
                                 control=list())
    if (status(modelEstimates) > 0){ 
        fake = rep(1000,length(d$tv))
        names(fake) <- names(d$tv)
        out <- fake #Livehack for ezsim to exclude failure results later
    }else{
        coef <- coefficients(modelEstimates)
        out <- c(coef$beta,coef$rhoY, coef$sigmaV, coef$sigmaU, coef$rhoV, coef$rhoU, coef$mu)
        
    }
    return(out)
}


#' @title True value for simulation
#'
#' @description
#' \code{spfrontier.true.value} returns true parameter values for a simulation process
#' 
#' @details
#' The \code{spfrontier.true.value} function should notbe used directly, it is exported for supporting \code{\link{ezsim}}
#' 
#' @rdname simulation

spfrontier.true.value <- function(){
    tv <- c(beta0, beta1, beta2)
    tvNames <- c("Beta0","Beta1","Beta2")
    if(!is.null(rhoY)){
        tv <- c(tv, rhoY)
        tvNames <- c(tvNames, "rhoY")
    }
    tv <- c(tv, sigmaV, sigmaU)
    tvNames <- c(tvNames, "SigmaV","SigmaU")
    if(!is.null(rhoV)){
        tv <- c(tv, rhoV)
        tvNames <- c(tvNames, "rhoV")
    }
    if(!is.null(rhoU)){
        tv <- c(tv, rhoU)
        tvNames <- c(tvNames, "rhoU")
    }
    if(!is.null(mu)){
        tv <- c(tv, mu)
        tvNames <- c(tvNames, "mu")
    }
    names(tv) <- tvNames
    return(tv)
}




#' @title Spatial stochastic frontier model simulation tests
#'
#' @description
#' \code{ezsimspfrontier} tests estimators of a spatial stochastic frontier model with different parameters
#' 
#' @details
#' The \code{ezsimspfrontier} function executes multiple calls of the \code{spfrontier} estimator on a simulated data set, 
#' generated on the base of provided parameters. The resulting estimates can be analysed for biasedness, efficiency, etc.
#' 
#' 
#' @param runs a number of simulated samples 
#' @param params a set with parameters to be used in simulation. There are predefined parameter sets:\cr
#' params000 - a non-spatial stochastic frontier with half-normal inefficiencies\cr
#' params000T - a non-spatial stochastic frontier with truncated normal inefficiencies\cr
#' params100 - a stochastic frontier with spatial lags of a dependent variable and with half-normal inefficiencies\cr
#' params100T - a stochastic frontier with spatial lags of a dependent variable and with truncated normal inefficiencies\cr
#' params110 -  stochastic frontier with spatial lags of a dependent variable and of a symmetric error component and with half-normal inefficiencies\cr
#' params111 -  stochastic frontier with spatial lags of a dependent variable, a symmetric error component, and an inefficiency error component and with half-normal inefficiencies\cr
#' params011 -  stochastic frontier with spatial lags of a symmetric error component and an inefficiency error component and with half-normal inefficiencies\cr
#' params001 -  stochastic frontier with spatial lags  an inefficiency error component and with half-normal inefficiencies\cr
#'     
#' @param autoSave save intermediate results to files. See \code{\link{ezsim}} for details.
#' @param seed a state for random number generation in R. If NULL (default), the initial state is random. See \code{\link{set.seed}} for details.
#' @param inefficiency sets the distribution for inefficiency error component. Possible values are 'half-normal' (for half-normal distribution) and 'truncated' (for truncated normal distribution). 
#' By default set to 'half-normal'. See references for explanations
#' @param logging an optional level of logging. Possible values are 'quiet','warn','info','debug'. 
#' By default set to quiet.
#' 
#' 
#' @keywords spatial stochastic frontier, simulation
#' @export
#' @seealso 
#' \code{\link{ezsim}}
#' 
#' @examples
#' 
#' # Define parameter values
#' # params <- list(n=c(50,100), 
#' #                  sigmaX=10, 
#' #                  beta0=1, 
#' #                  beta1=-2, 
#' #                  beta2=3, 
#' #                  sigmaV=0.2, 
#' #                  sigmaU=0.75)
#' 
#' # Run simulations (10 runs)
#' # res <- ezsimspfrontier(10, 
#' #                          params = params,  
#' #                          seed = 99, 
#' #                          inefficiency = "half-normal",
#' #                          logging = "quiet")
#' # 
#' # Summary of simulation results
#' # summary(res)
#' 
#' # Plot estimates' convergence to true values and estimates' density
#' #plot(res)
#' #plot(res, 'density')
#' 
#' @rdname simulation

ezsimspfrontier <- function(runs, 
                            autoSave = 0, 
                            params = list(n=c(50,100), sigmaX=10, beta0=1, beta1=-2, beta2=3, sigmaV=0.2, sigmaU=0.75),
                            seed = NULL,
                            inefficiency = "half-normal",
                            logging = "info"){
    if (!is.null(seed)) set.seed(seed)
    parDef <- createParDef(selection = params, banker = list(loggingLevel=logging,inefficiency=inefficiency,parDef=params))

    ezsim_spfrontier<-ezsim(
        m             = runs,
        run           = TRUE,
        parameter_def = parDef,
        dgp           = spfrontier.dgp,
        estimator     = spfrontier.estimator,
        true_value    = spfrontier.true.value,
        auto_save   = autoSave
    )
    
    ezsim_spfrontier <- clearFakes(ezsim_spfrontier)
    return(ezsim_spfrontier)
}


clearFakes = function(ezsim_ob){
    parSets = length(ezsim_ob$simulation_result)
    results = list()
    for (j in 1:parSets){
        results[[j]] = data.frame()
        runs = length(ezsim_ob$simulation_result[[j]])
        for (i in 1:runs){
            if (ezsim_ob$simulation_result[[j]][[runs-i+1]][1]==1000){
                ezsim_ob$simulation_result[[j]][[runs-i+1]] = NULL
            }else{
                results[[j]] = rbind(results[[j]],ezsim_ob$simulation_result[[j]][[runs-i+1]])
            }
        }
    }
    ezsim_ob = createSimulationTable(ezsim_ob)
    ezsim_ob$results = results
    return(ezsim_ob)
}