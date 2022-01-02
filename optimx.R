    out <- optimx(par,
                  gplvm.f,
                  gplvm.gr,
                  X=X,
                  groups=groups,
                  q=q,
                  KK=KK,
                  method="L-BFGS-B",
                  control=list(trace=T,
                               maximize=T,
                               kkt=FALSE,
                               maxit=1000,
                               starttests=FALSE)
    )
    
    
    
    
    
optimx<-function (par, fn, gr = NULL, hess = NULL, lower = -Inf, upper = Inf, 
method = c("Nelder-Mead", "BFGS"), itnmax = NULL, hessian = FALSE, 
control = list(), ...) 
{
    optcfg <- optimx.setup(par, fn, gr, hess, lower, upper, method, 
        itnmax, hessian, control, ...)
    if (optcfg$ctrl$starttests) {
        optchk <- optimx.check(par, optcfg$ufn, optcfg$ugr, optcfg$uhess, 
            lower, upper, hessian, optcfg$ctrl, have.bounds = optcfg$have.bounds, 
            usenumDeriv = optcfg$usenumDeriv, ...)
    }
    optcfg$ctrl$have.bounds <- optcfg$have.bounds
    if (!is.null(control$trace) && control$trace > 1) {
        cat("optcfg:")
        print(optcfg)
    }
    ansout <- optimx.run(par, optcfg$ufn, optcfg$ugr, optcfg$uhess, 
        lower, upper, optcfg$method, itnmax, hessian, optcfg$ctrl, 
        ...)
    details <- attr(ansout, "details")
    attr(ansout, "details") <- NULL
    if (optcfg$ctrl$maximize) {
        if (optcfg$ctrl$trace > 0) 
            cat("Reversing sign on objective, gradient, & hessian\n")
        ansout$value <- -ansout$value
        nlist <- dim(details)[[1]]
        for (i in 1:nlist) {
            details[[i, "ngatend"]] <- -details[[i, "ngatend"]]
            details[[i, "nhatend"]] <- -details[[i, "nhatend"]]
            details[[i, "hev"]] <- -details[[i, "hev"]]
        }
    }
    rownames(details) <- details[, "method"]
    ansout[, "kkt1"] <- as.logical(ansout[, "kkt1"])
    ansout[, "kkt2"] <- as.logical(ansout[, "kkt2"])
    answer <- structure(ansout, details = details, maximize = optcfg$ctrl$maximize, 
        npar = optcfg$npar, follow.on = optcfg$ctrl$follow.on, 
        class = c("optimx", "data.frame"))
    answer
}
    
    
    
    
    
        ansout <- optimx.run(par, fn, gr, hess=NULL, 
        lower, upper, method="L-BFGS-B", itnmax, hessian, optcfg$ctrl)
        
out<-optim(par,fn=gplvm.f,gr=gplvm.gr,X=X,groups=groups,q=q,KK=KK,control=list(maxit=1))
    
    
    
    
    
    
    