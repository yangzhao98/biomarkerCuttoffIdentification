rm(list = ls())
library(clinfun)
library(parallel)
library(R2jags)
library(mlr)
library(BMS)
library(MCMCglmm)
library(tidyverse)
library(HDInterval)
library(TeachingDemos)
library(overlapping)
library(emdbook)
library(FNN)
library(hdrcde)
library(tidyverse)
library(BTDesign)
library(bridgesampling)

fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

## Sample size calculation for binary endpoint
### For Binary Endpoint
## With pooled variance
bEventPooled <- function(alpha = 0.05, beta = 0.20, delta = 0.05, theta = 0, 
                         r = 1, pE = 0.12, pC = 0.12) {
  nC <- ceiling((qnorm(alpha)+qnorm(beta))^2/(theta + delta)^2*(pE*(1-pE)/r + pC*(1-pC)))
  nE <- nC*r
  n <- nC + nE
  return(list(nTotal = n, nTreatment = nE, nControl = nC))
}

bEventPooled(alpha = 0.025, beta = 0.2, delta = 0.10, theta = 0.05, 
             r = 2, pE = 0.70, pC = 0.75)

## Without pooled variance
bEventUnPooled <- function(alpha = 0.05, beta = 0.20, delta = 0.05, theta = 0, 
                           r = 2, pE = 0.12, pC = 0.12) {
  pbar <- (r*pE + pC)/(r + 1)
  nC <- ceiling((1+1/r)/(theta + delta)^2*(qnorm(alpha)*sqrt(pbar*(1-pbar)) + qnorm(beta)*sqrt((pE*(1-pE)/r+pC*(1-pC))/(1+1/r)))^2)
  nE <- nC*r
  n <- nC + nE
  return(list(nTotal = n, nTreatment = nE, nControl = nC))
}
bEventUnPooled(alpha = 0.025, beta = 0.2, delta = 0.05, theta = 0.0,
               r = 2, pE = 0.70, pC = 0.70)

## Simulation model
logitFunc <- function(a0 = -2.5, a1 = 0.1, a2 = 0.6, c = 3, 
                      p0 = 0.1, p1 = 0.3, bMin = 0, bMax = 5) {
  b <- seq(bMin, bMax, by = 0.05)
  logitp <- a0 + a1*b*(b < c) + a2*(b >= c)
  p <- exp(logitp)/(1 + exp(logitp))
  plot(x = b, y = p, col = "black", type = "l", 
       ylim = c(0,1.0), xlim = c(0,5), 
       xlab = "Biomarker", ylab = "p",
       main = bquote("logit(p) = " ~ .(a0) ~ "+" ~ .(a1) ~ "* B * I(B < " ~ 
                       .(c) ~ ") +" ~ .(a2) ~ "* B * I(B >= " ~ .(c) ~ ")"),
       cex.main = 0.9)
  abline(h = p0, col = "gray", lty = 2)
  abline(h = p1, col = "gray", lty = 2)
  abline(v = c, col = "red", lty = 2)
  points(x = b, y = p, pch = 1, col = "blue")
}
logitFunc(a0 = -3.0, a1 = 0.01, a2 = 2.0, c = 3, p0 = 0.05, p1 = 0.25)

### Identify optimal cutoff value for n-arms
SelectCutoffnArm <- function(nSim, n, bMin, bMax, 
                             a0, a1, a2, c, theta,
                             mean.Mu, var.Mu, tau.alpha, tau.beta,
                             ratioL, ratioH) {
  #print(nSim)
  set.seed(nSim)
  ## -- Section 1. Bayesian hierarchical model for basket trial 
  BHMNorm <- function(datList) {
    ## ---------------------------------------------
    ## -- dat
    ## --  (1) K = #. of baskets
    ## --  (2) y = #. of responses
    ## --  (3) n = #. of patients
    ## --  (4) zeta 
    Model <- "
    model{
    ## likelihood
    for (k in 1:(2*K)) {
    y[k] ~ dbin(p[k], n[k])
    logit(p[k]) <- rho[k]
    rho[k] ~ dnorm(mu[gIndex[k]], invtau[gIndex[k]])
    }
    ## prior
    for (g in 1:2){
    mu[g] ~ dnorm(mean.Mu, var.Mu)
    invtau[g] ~ dgamma(tau.alpha, tau.beta)
    }
    tau <- 1/invtau
    }
    "
    Fit <- jags(
      model.file = textConnection(Model),
      data = datList,
      parameters.to.save = c("p", "mu", "tau"),
      n.iter = 2000,
      n.burnin = 1000,
      n.chains = 1, n.thin = 1, progress.bar = "none")
    PoPp <- Fit$BUGSoutput$sims.list$p
    PoPMu<- Fit$BUGSoutput$sims.list$mu
    PoPTau <- Fit$BUGSoutput$sims.list$tau
    return(list(PoPp = PoPp, PoPMu = PoPMu, PoPTau = PoPTau))
  }
  
  ## Checking parameter
  K <- length(n)
  if (length(bMin) != K | length(bMax) != K | length(a0) != K | length(a1) != K |
      length(c) != K | length(theta) != K) {
    stop("Please, check the parameters!")
  }
  
  ## Generateing simulated dataset for biomarker cutoff identification
  datArray <- array(NA, dim = c(K,max(n),3))   ## Story the generated trial data
  Bio <- NULL
  for (k in 1:K) {
    B <- p <- y <- NULL
    #B <- runif(n[k], min = bMin[k], max = bMax[k])
    B <- rlnorm(n[k])
    p <- exp(a0[k] + a1[k]*B*(B <= c[k]) + a2[k]*(B > c[k]))/
      (1 + exp(a0[k] + a1[k]*B*(B <= c[k]) + a2[k]*(B > c[k])))
    y <- (runif(n[k]) <= p) + 0
    ## Story trial data
    datArray[k,1:n[k],1] <- B
    datArray[k,1:n[k],2] <- p
    datArray[k,1:n[k],3] <- y
    Bio <- c(Bio, B)
  }
  
  ##   
  quant <- seq(0.1, 0.9, by = 0.1)   ## Potential percentiles
  pCutoff <- matrix(0, nrow = 1, ncol = length(quant))
  colnames(pCutoff) <- paste(quant*100, "%", sep = "")
  pCutoff[,1:length(quant)] <- quantile(Bio, quant)    ## Potential cutoff values
  dat1Array <- array(NA, dim = c(K, length(quant), 5))  ## Story the observed data
  
  for (k in 1:K) {
    for (j in 1:length(quant)) {
      dat <- datArray[k,1:n[k],]
      dat1Array[k,j,1] <- sum(dat[dat[,1] <= pCutoff[j],3])    ## y0 | bio < c
      dat1Array[k,j,2] <- length(dat[dat[,1] <= pCutoff[j],3]) ## n0 | bio < c
      dat1Array[k,j,3] <- sum(dat[dat[,1] > pCutoff[j],3])   ## y1 | bio >= c
      dat1Array[k,j,4] <- length(dat[dat[,1] > pCutoff[j],3])## n1 | bio >= c
      dat1Array[k,j,5] <- length(dat[dat[,1] > pCutoff[j],3])/(length(dat[dat[,1] > pCutoff[j],3]) + length(dat[dat[,1] <= pCutoff[j],3]))
    }
  }
  
  ## 
  PoPMat <- matrix(NA, ncol = K*5, nrow = length(quant))## Story the posterior prob
  PoPMu <- PoPTau <- NULL
  row.names(PoPMat) <- paste(quant*100,"%",sep = "")
  colNames <- NULL
  for (k in 1:K) {
    colNames <- c(colNames, paste("p",k,"-", sep = ""), 
                  paste("p",k,"+",sep = ""),
                  paste("Pr",k,"(p+-p->",theta[k],")",sep = ""),
                  paste("n+",k,"/n-",k, sep = ""),
                  paste("l+",k,"/l-",k, sep = ""))
  }
  colnames(PoPMat) <- colNames
  for (j in 1:length(quant)) {
    ymat <- nmat <- NULL
    for (k in 1:K) {
      ymat <- c(ymat, c(dat1Array[k,j,1], dat1Array[k,j,3]))
      nmat <- c(nmat, c(dat1Array[k,j,2], dat1Array[k,j,4]))
    }
    #ymat <- c(dat1Array[,j,1],dat1Array[,j,3])
    #nmat <- c(dat1Array[,j,2],dat1Array[,j,4])
    PoP2 <- BHMNorm(datList = list(
      K = K,
      y = ymat,
      n = nmat,
      gIndex = rep(c(1,2),K),
      mean.Mu = mean.Mu, var.Mu = var.Mu,
      tau.alpha = tau.alpha, tau.beta = tau.beta))
    for (k in 1:K) {
      PoPMat[j,(5*k - 4):(5*k - 3)] <- colMeans(PoP2$PoPp[,c(2*k-1,2*k)])
      PoPMat[j,5*k - 2] <- mean((PoP2$PoPp[,k*2] - PoP2$PoPp[,k]) >= theta[k])
      PoPMat[j,5*k - 1] <- dat1Array[k,j,5]
      PoPMat[j,5*k] <- abs(log(mean(PoP2$PoPp[,2*k])/(1 - mean(PoP2$PoPp[,2*k]))) - log(theta[k]/(1 - theta[k])))/
        abs(log(mean(PoP2$PoPp[,2*k-1])/(1 - mean(PoP2$PoPp[,2*k-1]))) - log(theta[k]/(1 - theta[k])))
      #PoPMat[j,3*k] <- mean(PoP2$PoPp[,k+3] - PoP2$PoPp[,k] >= theta[k])
    }
    PoPMu <- c(PoPMu, mean(PoP2$PoPMu))
    PoPTau<- c(PoPTau, mean(PoP2$PoPTau))
  }
  
  ## Identify the optimal cutoff for each group
  index1Per <- index2Per <- index1 <- index2 <- NULL
  for (k in 1:K) {
    #ind <- (PoPMat[,4*k-2] - PoPMat[,4*k-3]) > theta[k]
    indPoP <- PoPMat[,4*k-1]
    ind2 <- PoPMat[,4*k] >= ratioL & PoPMat[,4*k] <= ratioH
    #ind[!ind2] <- FALSE
    indPoP[!ind2] <- 0
    #if (sum(ind) > 0) {
    #index2Per <- c(index2Per, quant[which.max(ind)]*10)
    index2Per <- c(index2Per, quant[which.max(indPoP)]*10)
    #index2 <- c(index2, pCutoff[which.max((PoPMat[,4*k-2] - PoPMat[,4*k-3]) > theta[k])]) 
    index2 <- c(index2, pCutoff[which.max(indPoP)])
    #} else {
    #  index2Per <- c(index2Per, quant[which.max(PoPMat[,4*k-3])]*10)
    #  index2 <- c(index2, pCutoff[which.max(PoPMat[,4*k-3])])      
    #}
    #index2Per <- c(index2Per, quant[which.max((PoPMat[,3*k-1]) > theta[k])]*10)
    #index2 <- c(index2, pCutoff[which.max((PoPMat[,3*k-1]) > theta[k])])
  }
  
  ## Summarize the results
  res <- NULL
  #res$Cutoff1 <- index1[which.max(index1)]
  #res$Index1 <- index1Per[which.max(index1)]
  #dat1Stage1 <- dat1Array[,index1Per[which.max(index1)],]
  #row.names(dat1Stage1) <- paste("Arm", 1:K, sep = "")
  #colnames(dat1Stage1) <- c("y1-","n1-","y1+","n1+")
  #dat1Stage1 <- as.data.frame(dat1Stage1)
  #res$dat1Stage1 <- dat1Stage1
  
  res$Cutoff2 <- index2[which.max(index2)]
  res$Index2 <- index2Per[which.max(index2)]
  dat2Stage1 <- dat1Array[,index2Per[which.max(index2)],] 
  row.names(dat2Stage1) <- paste("Arm", 1:K, sep = "")
  colnames(dat2Stage1) <- c("y1-","n1-","y1+","n1+","n+/n-")
  dat2Stage1 <- as.data.frame(dat2Stage1)
  res$dat2Stage1 <- dat2Stage1
  return(list(Result1 = res, Result2 = PoPMat, 
              PoPMu = PoPMu, PoPTau = PoPTau,
              pCutoff = pCutoff))
  }


SummarySelectCutoffnArm <- function(fit, nSimulation, cp = 6) {
  K <- nrow(fit[[1]]$Result1$dat2Stage1)
  nCP <- nrow(fit[[1]]$Result2)
  ResultMat <- matrix(0, nrow = nCP, ncol = K*5)
  ResultSel <- NULL
  PoPMu <- PoPTau <- NULL
  Cutoff1 <- Cutoff2 <- 0
  CutoffMat <- rep(0,nCP)
  #SelPercent <- NULL
  Stage1Dat <- Stage2Dat <- matrix(0, nrow = K, ncol = 5)
  for (i in 1:nSimulation) {
    ResultMat <- ResultMat + fit[[i]]$Result2/nSimulation
    ResultSel <- c(ResultSel, fit[[i]]$Result1$Index2)
    PoPMu <- c(PoPMu, fit[[i]]$PoPMu[cp])
    PoPTau<- c(PoPTau, fit[[i]]$PoPTau[cp])
    #Cutoff1 <- Cutoff1 + fit[[i]]$Result1$Cutoff1/nSimulation
    Cutoff2 <- Cutoff2 + fit[[i]]$Result1$Cutoff2/nSimulation
    CutoffMat <- CutoffMat + fit[[i]]$pCutoff/nSimulation
    #cutpoint <- pmax(which.max(fit[[i]]$Result2[,3]),
    #                 which.max(fit[[i]]$Result2[,3*2]),
    #                 which.max(fit[[i]]$Result2[,3*3]))
    #SelPercent <- c(SelPercent, cutpoint)
    #Stage1Dat <- Stage1Dat + fit[[i]]$Result1$dat1Stage1/nSimulation
    Stage2Dat <- Stage2Dat + fit[[i]]$Result1$dat2Stage1/nSimulation
  }
  names(Cutoff1) <- NULL
  names(Cutoff2) <- NULL
  Stage1Dat <- as.data.frame(Stage1Dat)
  Stage2Dat <- as.data.frame(Stage2Dat)
  ResultMat <- as.data.frame(ResultMat)
  SelMat <- table(ResultSel)/nSimulation
  return(list(ResultMat = ResultMat, SelMat = SelMat,
              PoP1Mu = PoPMu, PoP1Tau = PoPTau,
              #Cutoff1 = Cutoff1, 
              Cutoff2 = Cutoff2,
              CutoffMat = CutoffMat, 
              #Stage1Dat = Stage1Dat,
              Stage2Dat = Stage2Dat))
}

n <- c(30, 30, 30)
bMin <- c(0,0,0)
bMax <- c(5,5,5)
theta <- c(0.15,0.15,0.15)#c(0.3,0.3,0.3) + 0.1
nSimulation <- 500
tau.alpha <- tau.beta <- 0.001
ratioL <- 0.2
ratioH <- 0.8

## -- Scenario 1. 
logitFunc(a0 = -3.0, a1 = 0.01, a2 = 2, c = 0.65, p0 = 0.05, p1 = 0.25)
## - Stage 1. Identify the optimal cutoff value
p0 <- c(0.05,0.05,0.05)
p1 <- c(0.25,0.25,0.25)
a0 <- c(-3.0,-3.0,-3.0)
a1 <- c(0.01,0.01,0.01)
a2 <- c(0.65,0.65,0.65)
c <- c(3.0,3.0,3.0)
SimCutoffnArm1 <- mclapply(1:nSimulation,
                           function(nSim) SelectCutoffnArm(nSim = nSim, n = n, 
                                                           bMin = bMin, bMax = bMax,
                                                           a0 = a0, a1 = a1, a2 = a2, 
                                                           c = c, theta = theta,
                                                           mean.Mu = log(0.05/0.95), var.Mu = 10,
                                                           tau.alpha = 0.001, tau.beta = 0.001,
                                                           ratioL = ratioL, ratioH = ratioH),
                           mc.cores = 4)
SimCPRes1 <- SummarySelectCutoffnArm(fit = SimCutoffnArm1, 
                                     nSimulation = nSimulation, cp = c[1]/5*10)



##  Bayesian Model Selection Framework
SelectCutoffnArmBMS <- function(nSim, n, bMin, bMax, 
                               a0, a1, a2, c, theta,
                               mean.Mu, var.Mu, tau.alpha, tau.beta,
                               ratioL, ratioH) {
  set.seed(nSim)
  ## -- Marginal likelihood function of Mj, L(D|Mj)
  # marlikM <- function(n,y) {
  #   tauInt <- Vectorize(function(tau){
  #     muInt <- Vectorize(function(mu){
  #       pInt <- Vectorize(function(p){
  #         sum(dbinom(y,n,p, log = TRUE)*dnorm(log(p/(1-p)),mu,tau), log = TRUE)*
  #           dnorm(mu,log(0.3/0.7),10, log = TRUE)*dgamma(tau, 2, 20, log = TRUE)
  #       })
  #       integrate(pInt, lower = 0, upper = 1, subdivisions = 2000)$value
  #     })
  #     integrate(muInt, lower = -Inf, upper = Inf, subdivisions = 2000)$value
  #   })
  #   ml <- integrate(tauInt, lower = 0, upper = Inf, subdivisions = 2000)$value
  #   return(ml)
  # }
  # marlikM(n = nmat, y = ymat)
  
  ## -------
  ## Step 1.
  ## Bayesian hierarchical model
  getSampleModel <- function(data, niter = 2000, nburn = 1000) {
    model <- "
    model{
      for (k in 1:(2*K)) {
        y[k] ~ dbin(p[k], n[k])
        logit(p[k]) <- rho[k]
        rho[k] ~ dnorm(mu[gIndex[k]], invtau[gIndex[k]])
      }
      for (g in 1:2) {
        mu[g] ~ dnorm(mean.Mu, var.Mu)
        invtau[g] ~ dgamma(tau.alpha, tau.beta)
      }
      tau <- 1/invtau
    }
    "
    fit <- jags(
      model.file = textConnection(model),
      data = data,
      parameters.to.save = c("p", "mu", "tau"),
      n.iter = niter,
      n.burnin = nburn,
      n.chains = 1, n.thin = 1, progress.bar = "none")
    return(fit)
  }
  
  ## Step 2. Specify the unnormalized log posterior function
  logPosteriorModel <- function(samples.row, data) {
    mu1 <- samples.row[paste0("mu[",unique(data$gIndex)[1],"]", sep = "")]
    mu2 <- samples.row[paste0("mu[",unique(data$gIndex)[2],"]", sep = "")]
    invtau1 <- samples.row[paste0("tau[",unique(data$gIndex)[1],"]", sep = "")]
    invtau2 <- samples.row[paste0("tau[",unique(data$gIndex)[2],"]", sep = "")]
    p1 <- samples.row[paste0("p[",seq(data$K)*2-1, "]", sep = "")]
    p2 <- samples.row[paste0("p[",seq(data$K)*2, "]", sep = "")] 
    y1 <- data$y[seq(data$K)*2-1]
    y2 <- data$y[seq(data$K)*2]
    n1 <- data$n[seq(data$K)*2-1]
    n2 <- data$n[seq(data$K)*2]
    
    sum(dbinom(y1, n1, p1, log = TRUE)) +
      sum(dnorm(log(p1/(1-p1)), mu1, invtau1, log = TRUE)) +
      sum(dbinom(y2, n2, p2, log = TRUE)) +
      sum(dnorm(log(p2/(1-p2)), mu2, invtau2, log = TRUE)) +
      dnorm(mu1, data$mean.Mu, data$var.Mu, log = TRUE) +
      dnorm(mu2, data$mean.Mu, data$var.Mu, log = TRUE) +
      dgamma(invtau1, data$tau.alpha, data$tau.beta, log = TRUE) +
      dgamma(invtau2, data$tau.alpha, data$tau.beta, log = TRUE)
  }
  
  # ## Step 3. Specify the paarameter bounds
  # cn <- colnames(sampleModel$BUGSoutput$sims.matrix)
  # cn <- cn[cn != "deviance"]
  # lb <- rep(-Inf, length(cn))
  # ub <- rep(Inf, length(cn))
  # names(lb) <- names(ub) <- cn
  # lb[c("tau[1]","tau[2]")] <- 0
 
  ## Step 4. COmpute the marginal likelihoods
  # Modelbridge <- bridge_sampler(samples = sampleModel,
  #                               data = datList,
  #                               log_posterior = logPosteriorModel,
  #                               lb = lb, ub = ub, silent = TRUE)
  
  ## Checking parameter
  K <- length(n)
  if (length(bMin) != K | length(bMax) != K | length(a0) != K | length(a1) != K |
      length(c) != K | length(theta) != K) {
    stop("Please, check the parameters!")
  }
  ## Generateing simulated dataset for biomarker cutoff identification
  datArray <- array(NA, dim = c(K,max(n),3))   ## Story the generated trial data
  Bio <- NULL
  for (k in 1:K) {
    B <- p <- y <- NULL
    B <- runif(n[k], min = bMin[k], max = bMax[k])
    #B <- rlnorm(n[k])
    p <- exp(a0[k] + a1[k]*B*(B <= c[k]) + a2[k]*(B > c[k]))/
      (1 + exp(a0[k] + a1[k]*B*(B <= c[k]) + a2[k]*(B > c[k])))
    y <- (runif(n[k]) <= p) + 0
    ## Story trial data
    datArray[k,1:n[k],1] <- B
    datArray[k,1:n[k],2] <- p
    datArray[k,1:n[k],3] <- y
    Bio <- c(Bio, B)
  }
  
  ##   
  quant <- seq(0.1, 0.9, by = 0.1)   ## Potential percentiles
  pCutoff <- matrix(0, nrow = 1, ncol = length(quant))
  colnames(pCutoff) <- paste(quant*100, "%", sep = "")
  pCutoff[,1:length(quant)] <- quantile(Bio, quant)    ## Potential cutoff values
  dat1Array <- array(NA, dim = c(K, length(quant), 5))  ## Story the observed data
  
  for (k in 1:K) {
    for (j in 1:length(quant)) {
      dat <- datArray[k,1:n[k],]
      dat1Array[k,j,1] <- sum(dat[dat[,1] <= pCutoff[j],3])    ## y0 | bio < c
      dat1Array[k,j,2] <- length(dat[dat[,1] <= pCutoff[j],3]) ## n0 | bio < c
      dat1Array[k,j,3] <- sum(dat[dat[,1] > pCutoff[j],3])   ## y1 | bio >= c
      dat1Array[k,j,4] <- length(dat[dat[,1] > pCutoff[j],3])## n1 | bio >= c
      dat1Array[k,j,5] <- length(dat[dat[,1] > pCutoff[j],3])/(length(dat[dat[,1] > pCutoff[j],3]) + length(dat[dat[,1] <= pCutoff[j],3]))
    }
  }
  
  ## 
  PoPMat <- matrix(NA, ncol = K*5, nrow = length(quant))## Story the posterior prob
  PoPMu <- PoPTau <- NULL
  row.names(PoPMat) <- paste(quant*100,"%",sep = "")
  colNames <- NULL
  for (k in 1:K) {
    colNames <- c(colNames, paste("p",k,"-", sep = ""), 
                  paste("p",k,"+",sep = ""),
                  paste("Pr",k,"(p+-p->",theta[k],")",sep = ""),
                  paste("n+",k,"/n-",k, sep = ""),
                  paste("l+",k,"/l-",k, sep = ""))
  }
  colnames(PoPMat) <- colNames
  pMD <- NULL
  for (j in 1:length(quant)) {
    ymat <- nmat <- NULL
    for (k in 1:K) {
      ymat <- c(ymat, c(dat1Array[k,j,1], dat1Array[k,j,3]))
      nmat <- c(nmat, c(dat1Array[k,j,2], dat1Array[k,j,4]))
    }
    datList = list(
      K = K,
      y = ymat,
      n = nmat,
      gIndex = rep(c(1,2),K),
      mean.Mu = 0, var.Mu = 10,
      tau.alpha = 0.01, tau.beta = 0.01)
    sampleModel <- getSampleModel(datList)
    cn <- colnames(sampleModel$BUGSoutput$sims.matrix)
    cn <- cn[cn != "deviance"]
    lb <- rep(-Inf, length(cn))
    ub <- rep(Inf, length(cn))
    names(lb) <- names(ub) <- cn
    lb[c("tau[1]","tau[2]")] <- 0
    Modelbridge <- bridge_sampler(samples = sampleModel,
                                  data = datList,
                                  log_posterior = logPosteriorModel,
                                  lb = lb, ub = ub, silent = TRUE)
    pMD <- c(pMD, Modelbridge$logml)
  }
  
  ## Identify the optimal cutoff for each group
  index1Per <- index2Per <- index1 <- index2 <- NULL
  for (k in 1:K) {
    ## Sample size ratio
    indPoP <- PoPMat[,5*k-1]
    ind2 <- PoPMat[,4*k] >= ratioL & PoPMat[,4*k] <= ratioH
    indPoP[!ind2] <- 0
    index2Per <- c(index2Per, quant[which.max(indPoP)]*10)
    index2 <- c(index2, pCutoff[which.max(indPoP)])
    ## First meet (p+ > p- + delta)
    ind1 <- ((PoPMat[,5*k-3] - PoPMat[,5*k-4]) > theta[k])
    index1Per <- c(index1Per, quant[which.max(ind1)]*10)
    index1 <- c(index1, pCutoff[which.max(ind1)])
  }
  
  ## Summarize the results
  res <- NULL
  #res$Cutoff1 <- index1[which.max(index1)]
  #res$Index1 <- index1Per[which.max(index1)]
  #dat1Stage1 <- dat1Array[,index1Per[which.max(index1)],]
  #row.names(dat1Stage1) <- paste("Arm", 1:K, sep = "")
  #colnames(dat1Stage1) <- c("y1-","n1-","y1+","n1+")
  #dat1Stage1 <- as.data.frame(dat1Stage1)
  #res$dat1Stage1 <- dat1Stage1
  
  res$Cutoff2 <- index2[which.max(index2)]
  res$Index2 <- index2Per[which.max(index2)]
  dat2Stage1 <- dat1Array[,index2Per[which.max(index2)],] 
  row.names(dat2Stage1) <- paste("Arm", 1:K, sep = "")
  colnames(dat2Stage1) <- c("y1-","n1-","y1+","n1+","n+/n-")
  dat2Stage1 <- as.data.frame(dat2Stage1)
  res$dat2Stage1 <- dat2Stage1
  return(list(Result1 = res, Result2 = PoPMat, 
              PoPMu = PoPMu, PoPTau = PoPTau,
              pCutoff = pCutoff))
  }


SelectCutoffnArmBMSLogistic <- function(nSim, n, bMin, bMax, 
                                a0, a1, a2, c, theta,
                                mean.Mu, var.Mu, tau.alpha, tau.beta,
                                ratioL, ratioH) {
  set.seed(nSim)

  ## -------
  ## Step 1.
  ## Bayesian hierarchical model
  getSampleModel <- function(data, niter = 2000, nburn = 1000) {
    model <- "
      model{
        for (i in 1:N) {
          ind[i] <- ifelse(B[i] <= C, 1, 0)
          y[i] ~ dbern(p[i])
          logit(p[i]) <- beta0[arm[i]] + beta1[arm[i]]*B[i]*ind[i] 
                                       + beta2[arm[i]]*B[i]*ind[i] 
        }
      
        for (k in 1:K) {
          beta0[k] ~ dnorm(mu, 10)
          beta1[k] ~ dnorm(0, 10)
          beta2[k] ~ dnorm(0, 10)
        }
        mu ~ dnorm(0, 10)
      }
    "    
    fit <- jags(
      model.file = textConnection(model),
      data = data,
      parameters.to.save = c("beta0", "beta1", "beta2"),
      n.iter = niter,
      n.burnin = nburn,
      n.chains = 1, n.thin = 1, progress.bar = "none")
    return(fit)
  }
  
  ## Step 2. Specify the unnormalized log posterior function
  logPosteriorModel <- function(samples.row, data) {
    mu1 <- samples.row[paste0("mu[",unique(data$gIndex)[1],"]", sep = "")]
    mu2 <- samples.row[paste0("mu[",unique(data$gIndex)[2],"]", sep = "")]
    invtau1 <- samples.row[paste0("tau[",unique(data$gIndex)[1],"]", sep = "")]
    invtau2 <- samples.row[paste0("tau[",unique(data$gIndex)[2],"]", sep = "")]
    p1 <- samples.row[paste0("p[",seq(data$K)*2-1, "]", sep = "")]
    p2 <- samples.row[paste0("p[",seq(data$K)*2, "]", sep = "")] 
    y1 <- data$y[seq(data$K)*2-1]
    y2 <- data$y[seq(data$K)*2]
    n1 <- data$n[seq(data$K)*2-1]
    n2 <- data$n[seq(data$K)*2]
    
    sum(dbinom(y1, n1, p1, log = TRUE)) +
      sum(dnorm(log(p1/(1-p1)), mu1, invtau1, log = TRUE)) +
      sum(dbinom(y2, n2, p2, log = TRUE)) +
      sum(dnorm(log(p2/(1-p2)), mu2, invtau2, log = TRUE)) +
      dnorm(mu1, data$mean.Mu, data$var.Mu, log = TRUE) +
      dnorm(mu2, data$mean.Mu, data$var.Mu, log = TRUE) +
      dgamma(invtau1, data$tau.alpha, data$tau.beta, log = TRUE) +
      dgamma(invtau2, data$tau.alpha, data$tau.beta, log = TRUE)
  }
  
  # ## Step 3. Specify the paarameter bounds
  # cn <- colnames(sampleModel$BUGSoutput$sims.matrix)
  # cn <- cn[cn != "deviance"]
  # lb <- rep(-Inf, length(cn))
  # ub <- rep(Inf, length(cn))
  # names(lb) <- names(ub) <- cn
  # lb[c("tau[1]","tau[2]")] <- 0
  
  ## Step 4. COmpute the marginal likelihoods
  # Modelbridge <- bridge_sampler(samples = sampleModel,
  #                               data = datList,
  #                               log_posterior = logPosteriorModel,
  #                               lb = lb, ub = ub, silent = TRUE)
  
  ## Checking parameter
  K <- length(n)
  if (length(bMin) != K | length(bMax) != K | length(a0) != K | length(a1) != K |
      length(c) != K | length(theta) != K) {
    stop("Please, check the parameters!")
  }
  ## Generateing simulated dataset for biomarker cutoff identification
  datArray <- array(NA, dim = c(K,max(n),3))   ## Story the generated trial data
  Bio <- NULL
  for (k in 1:K) {
    B <- p <- y <- NULL
    B <- runif(n[k], min = bMin[k], max = bMax[k])
    p <- exp(a0[k] + a1[k]*B*(B <= c[k]) + a2[k]*B*(B > c[k]))/
      (1 + exp(a0[k] + a1[k]*B*(B <= c[k]) + a2[k]*B*(B > c[k])))
    y <- (runif(n[k]) <= p) + 0
    ## Story trial data
    datArray[k,1:n[k],1] <- B
    datArray[k,1:n[k],2] <- p
    datArray[k,1:n[k],3] <- y
    Bio <- c(Bio, B)
  }
  
  ##   
  quant <- seq(0.1, 0.9, by = 0.1)   ## Potential percentiles
  pCutoff <- matrix(0, nrow = 1, ncol = length(quant))
  colnames(pCutoff) <- paste(quant*100, "%", sep = "")
  pCutoff[,1:length(quant)] <- quantile(Bio, quant)    ## Potential cutoff values
  dat1Array <- array(NA, dim = c(K, length(quant), 5))  ## Story the observed data
  
  dat2 <- NULL
  for (k in 1:K) {
    #for (j in 1:length(quant)) {
      dat <- as.data.frame(datArray[k,1:n[k],])
      dat1 <- cbind(dat, arm = k)
      dat2 <- rbind(dat2, dat1)
      #dat1Array[k,j,1] <- sum(dat[dat[,1] <= pCutoff[j],3])    ## y0 | bio < c
      #dat1Array[k,j,2] <- length(dat[dat[,1] <= pCutoff[j],3]) ## n0 | bio < c
      #dat1Array[k,j,3] <- sum(dat[dat[,1] > pCutoff[j],3])     ## y1 | bio >= c
      #dat1Array[k,j,4] <- length(dat[dat[,1] > pCutoff[j],3])  ## n1 | bio >= c
      #dat1Array[k,j,5] <- length(dat[dat[,1] > pCutoff[j],3])/(length(dat[dat[,1] > pCutoff[j],3]) + length(dat[dat[,1] <= pCutoff[j],3]))
    #}
  }
  datList <- list(N = length(dat2$arm),
                  K = max(unique(dat2$arm)),
                  y = dat2$V3,
                  B = dat2$V1,
                  arm = dat2$arm, C = 3)
  
  sampleModel <- getSampleModel(data = datList)
  
  for (k in 1:K) {
    for (j in 1:length(quant)) {
      dat <- datArray[k,1:n[k],]
      dat1Array[k,j,1] <- sum(dat[dat[,1] <= pCutoff[j],3])    ## y0 | bio < c
      dat1Array[k,j,2] <- length(dat[dat[,1] <= pCutoff[j],3]) ## n0 | bio < c
      dat1Array[k,j,3] <- sum(dat[dat[,1] > pCutoff[j],3])   ## y1 | bio >= c
      dat1Array[k,j,4] <- length(dat[dat[,1] > pCutoff[j],3])## n1 | bio >= c
      dat1Array[k,j,5] <- length(dat[dat[,1] > pCutoff[j],3])/(length(dat[dat[,1] > pCutoff[j],3]) + length(dat[dat[,1] <= pCutoff[j],3]))
    }
  }
  
  ## 
  PoPMat <- matrix(NA, ncol = K*5, nrow = length(quant))## Story the posterior prob
  PoPMu <- PoPTau <- NULL
  row.names(PoPMat) <- paste(quant*100,"%",sep = "")
  colNames <- NULL
  for (k in 1:K) {
    colNames <- c(colNames, paste("p",k,"-", sep = ""), 
                  paste("p",k,"+",sep = ""),
                  paste("Pr",k,"(p+-p->",theta[k],")",sep = ""),
                  paste("n+",k,"/n-",k, sep = ""),
                  paste("l+",k,"/l-",k, sep = ""))
  }
  colnames(PoPMat) <- colNames
  pMD <- NULL
  for (j in 1:length(quant)) {
    ymat <- nmat <- NULL
    for (k in 1:K) {
      ymat <- c(ymat, c(dat1Array[k,j,1], dat1Array[k,j,3]))
      nmat <- c(nmat, c(dat1Array[k,j,2], dat1Array[k,j,4]))
    }
    datList = list(
      K = K,
      y = ymat,
      n = nmat,
      gIndex = rep(c(1,2),K),
      mean.Mu = 0, var.Mu = 10,
      tau.alpha = 2, tau.beta = 20)
    sampleModel <- getSampleModel(datList)
    cn <- colnames(sampleModel$BUGSoutput$sims.matrix)
    cn <- cn[cn != "deviance"]
    lb <- rep(-Inf, length(cn))
    ub <- rep(Inf, length(cn))
    names(lb) <- names(ub) <- cn
    lb[c("tau[1]","tau[2]")] <- 0
    Modelbridge <- bridge_sampler(samples = sampleModel,
                                  data = datList,
                                  log_posterior = logPosteriorModel,
                                  lb = lb, ub = ub, silent = TRUE)
    pMD <- c(pMD, Modelbridge$logml)
  }
  
  ## Identify the optimal cutoff for each group
  index1Per <- index2Per <- index1 <- index2 <- NULL
  for (k in 1:K) {
    #ind <- (PoPMat[,4*k-2] - PoPMat[,4*k-3]) > theta[k]
    indPoP <- PoPMat[,4*k-1]
    ind2 <- PoPMat[,4*k] >= ratioL & PoPMat[,4*k] <= ratioH
    #ind[!ind2] <- FALSE
    indPoP[!ind2] <- 0
    #if (sum(ind) > 0) {
    #index2Per <- c(index2Per, quant[which.max(ind)]*10)
    index2Per <- c(index2Per, quant[which.max(indPoP)]*10)
    #index2 <- c(index2, pCutoff[which.max((PoPMat[,4*k-2] - PoPMat[,4*k-3]) > theta[k])]) 
    index2 <- c(index2, pCutoff[which.max(indPoP)])
    #} else {
    #  index2Per <- c(index2Per, quant[which.max(PoPMat[,4*k-3])]*10)
    #  index2 <- c(index2, pCutoff[which.max(PoPMat[,4*k-3])])      
    #}
    #index2Per <- c(index2Per, quant[which.max((PoPMat[,3*k-1]) > theta[k])]*10)
    #index2 <- c(index2, pCutoff[which.max((PoPMat[,3*k-1]) > theta[k])])
  }
  
  ## Summarize the results
  res <- NULL
  #res$Cutoff1 <- index1[which.max(index1)]
  #res$Index1 <- index1Per[which.max(index1)]
  #dat1Stage1 <- dat1Array[,index1Per[which.max(index1)],]
  #row.names(dat1Stage1) <- paste("Arm", 1:K, sep = "")
  #colnames(dat1Stage1) <- c("y1-","n1-","y1+","n1+")
  #dat1Stage1 <- as.data.frame(dat1Stage1)
  #res$dat1Stage1 <- dat1Stage1
  
  res$Cutoff2 <- index2[which.max(index2)]
  res$Index2 <- index2Per[which.max(index2)]
  dat2Stage1 <- dat1Array[,index2Per[which.max(index2)],] 
  row.names(dat2Stage1) <- paste("Arm", 1:K, sep = "")
  colnames(dat2Stage1) <- c("y1-","n1-","y1+","n1+","n+/n-")
  dat2Stage1 <- as.data.frame(dat2Stage1)
  res$dat2Stage1 <- dat2Stage1
  return(list(Result1 = res, Result2 = PoPMat, 
              PoPMu = PoPMu, PoPTau = PoPTau,
              pCutoff = pCutoff))
  }
