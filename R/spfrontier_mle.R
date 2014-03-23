Infin <- -1e16
funcLogL <- function(parameters, quiet = F){
    
    y <- envirGet("y")
    X <- envirGet("X")
    k <- ncol(X)
    W_y <- envirGet("W_y")
    W_v <- envirGet("W_v")
    W_u <- envirGet("W_u")
    inefficiency <- envirGet("inefficiency")
    
    isSpY <- !is.null(W_y)
    isSpV <- !is.null(W_v)
    isSpU <- !is.null(W_u)
    isTN <- (inefficiency == "truncated")
          
    p <- paramsFromVector(parameters, k, isSpY, isSpV, isSpU, isTN)
    if(!quiet){
        pbeta <- olsenReparamBack(p)
        paramsToVector(pbeta, olsen = F)
        counter = envirCounter("funcLogL")
        logging(paste("Evaluating 'funcLogL' for parameters","[run=",counter,"]:"),paramsToVector(pbeta, olsen = F)) 
    }  
    n <- length(y)
    I <- diag(n)
    ret <- -Inf
    if ((p$nu>0) && (p$lambda>0) && (is.null(p$rhoY) || abs(p$rhoY)<1) 
        && (is.null(p$rhoV) || abs(p$rhoV)<1) 
        && (is.null(p$rhoU) || abs(p$rhoU)<1)){
        
        e <- y - X %*% p$gamma/p$nu
        
        if (isSpY){
            e <- e - p$rhoY * W_y %*% y
        }
        
        SpV <- I
        SpU <- I
        
        tryCatch({
            if (isSpV)
                SpV <- solve(I-p$rhoV*W_v)
            if (isSpU)
                SpU <- solve(I-p$rhoU*W_u)
            
            sigmaV <- 1/(p$nu*sqrt(1+p$lambda^2))
            sigmaU <- p$lambda * sigmaV
            
            mSigma = sigmaV^2*SpV%*%t(SpV)
            mOmega = sigmaU^2*SpU%*%t(SpU)
            mC = mSigma + mOmega
            imC <- solve(mC)
            mB <- mOmega%*%imC%*%mSigma 
            
            #Livehack for calculation precision
            rownames(mB) <- colnames(mB)
            mB <- as.matrix(nearPD(mB)$mat)
            
            mA <- mB %*% solve(mSigma)
            mD <- -mB %*% solve(mOmega)
            colnames(mB) <-rownames(mB)
            mu <- 0
            if(isTN){
                mu <- p$mu
            }
            vMu = rep(mu, n)
            r <- 0
            if (isSpV || isSpU){
                v <- pmvnorm(lower=rep(-Inf, n),upper = 0, mean=-vMu, sigma=mOmega)
                if (v<1e-200){
                    stop("Log-likelihood function value can not be calculated due to precision limits")
                }
                r <- -log(v)
                #logging("1:",r)
                r <- r + dmvnorm(x=as.vector(e),mean=vMu, sigma=mC,log=TRUE)
                
                #logging("2:",r)
                r <- r + log(pmvnorm(lower=rep(-Inf, n),upper = as.vector(mOmega%*%imC%*%(e-vMu)), mean=-vMu, sigma=mB))
                
                #logging("3:",r)
            }else{
                r <- r - n*log(pnorm(mu/sigmaU))
                r <- r + n * log(p$nu) 
                r <- r - 0.5 * p$nu^2*sum((e + vMu)^2)
                r <- r + sum(log(pnorm(mu*p$nu/p$lambda - p$lambda*e*p$nu)))
            }
            if (isSpY)
                r <- r + determinant(I-p$rhoY*W_y,logarithm = TRUE)$modulus
            ret <- r
        }, error = function(e){
            logging(e$message, level="warn")
        })
    }else{
        if(!quiet) logging("Parameters are out of space")
    }
    if(!quiet) logging(paste("funcLogL =",-ret,"\n"))
    if (ret == -Inf || is.nan(ret)) ret<- Infin # For fake finite differences
    return(-ret)
}

funcLogL.direct <- function(parameters){
    y <- envirGet("y")
    X <- envirGet("X")
    k <- ncol(X)
    W_y <- envirGet("W_y")
    W_v <- envirGet("W_v")
    W_u <- envirGet("W_u")
    inefficiency <- envirGet("inefficiency")
    
    isSpY <- !is.null(W_y)
    isSpV <- !is.null(W_v)
    isSpU <- !is.null(W_u)
    isTN <- (inefficiency == "truncated")
    
    p <- paramsFromVector(parameters, k, isSpY, isSpV, isSpU, isTN, olsen = F)
    polsen <- olsenReparam(p)
    return(funcLogL(paramsToVector(polsen), quiet = F))
}

funcGradient<-function(parameters){
    logging("Evaluating 'funcGradient' for parameters:",parameters)
    y = envirGet("y")
    X = envirGet("X")
    k <- ncol(X)
    W_y <- envirGet("W_y")
    W_v <- envirGet("W_v")
    W_u <- envirGet("W_u")
    inefficiency <- envirGet("inefficiency")
    
    isSpY <- !is.null(W_y)
    isSpV <- !is.null(W_v)
    isSpU <- !is.null(W_u)
    isTN <- (inefficiency == "truncated")
    
    p <- paramsFromVector(parameters, k, isSpY, isSpV, isSpU, isTN)
    
    yst <- y
    if(isSpY){
        yst <- y-p$rhoY*W_y%*%y
    }
    mu <- 0
    if (isTN){
        mu <- p$mu
    }
    omega <- p$nu * yst - X%*%p$gamma
    a <- -omega * p$lambda + mu*p$nu/p$lambda
    delta <- dnorm(a)/pnorm(a)
    N <- length(y)
    
    sigmaU <- p$lambda/(p$nu*sqrt(1+p$lambda^2))
    a2 <- mu/sigmaU
    delta2 <- (dnorm(a2)/pnorm(a2))
    
    dLdGamma <- t(omega+mu*p$nu)%*%X + t(delta)%*%X * p$lambda
    dLdNu <- -t(omega+mu*p$nu)%*%(y+mu) - t(delta)%*%(y*p$lambda-mu/p$lambda) + N * 1/p$nu - N*delta2*a2/p$nu
    dLdLambda <- -t(delta)%*%(omega+mu*p$nu/p$lambda^2)+N*delta2*mu*p$nu/(p$lambda^2*sqrt(1+p$lambda^2))
    grad = c(dLdGamma)
    names(grad) = paste("dGamma", seq(length(dLdGamma)), sep = "")
    if(isSpY){
        mat = solve(diag(N)-p$rhoY*W_y) %*% (-W_y)
        dLdRhoY = p$nu*(t(omega)+p$lambda*t(delta))%*%W_y%*%y+sum(diag(mat))
        ns <- names(grad)
        grad <- c(grad, dLdRhoY)
        names(grad) <- c(ns, "dRhoY")
    }
    ns <- names(grad)
    grad = c(grad, dLdNu,dLdLambda)
    names(grad) = c(ns, "dNu", "dLambda")
    if(isTN){
        dLdMu <- -N*delta2/sigmaU - p$nu*sum(omega+mu*p$nu)+sum(delta)*p$nu/p$lambda
        #dLdMu = -(mu - 1)
        ns <- names(grad)
        grad <- c(grad, dLdMu)
        names(grad) <- c(ns, "dMu")
    }
    return(-grad)
}