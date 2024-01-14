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

DKL <- function(alpha, beta, y1, n1, y0, n0) {
  Func <- function(x, y0, n0, y1, n1) {
    return(x^(y1 - y0)*(1 - x)^(n1 - y1 - n0 + y0))
  }
  KL <- log(beta(alpha + y0, beta + n0 - y0)/beta(alpha + y1, beta + n1 - y1)) +
    integrate(Func, lower = 0, upper = 1, y1 = y1, n1 = n1, y0 = y0, n0 = n0)$value
  return(KL)
}

PoPD <- function(alpha1, beta1, alpha0, beta0,
                 y1, n1, y0, n0) {
  intFunc <- function(y1, n1, y0, n0, pplus, pmin) {
    pplus^(alpha1 + y1 - 1)*(1 - pplus)^(beta1 + n1 - y1 - 1)*
      pmin^(alpha0 + y0 - 1)*(1 - pmin)^(beta0 + n0 - y0 - 1)/
      (betabn(x = y1, n = n1, a = alpha1, b = beta1)*
         betabn(x = y0, n = n0, a = alpha0, b = beta0))
  }

  intValue <- 1 - integrate(function(pmin) {
    sapply(pmin, function(pmin) {
      integrate(function(pplus)
        intFunc(y1 = y1, n1 = n1, y0 = y0, n0 = n0, pmin, pplus),
        lower = pmin, upper = 1)$value
    })
  }, lower = 0, upper = 1)$value

  return(intValue)
}

n <- 60
bMin <- 0; bMax <- 5
p0 <- 0.1
p1 <- 0.3
a0 <- -2
a1 <- 0.1
a2 <- 0.4
#a1 <-  log(p0/(1 - p0)) - a0
#a2 <-  log(p1/(1 - p1)) - a0
c <-  2.5


## Generating simulated data set
GenTrialData <- function(n, bMin, bMax,a0, a1, a2, c){
  B <- runif(n, min = bMin, max = bMax)
  c <- rep(c, n)
  p <- exp(a0 + a1*B*(B <= c) + a2*B*(B > c))/(1 + exp(a0 + a1*B*(B <= c) + a2*B*(B > c)))
  dat <- data.frame(
    p = p, biomarker = B
  )
  return(dat)
}
# dat <- GenTrialData(n = n, bMin = bMin, bMax = bMax, a0 = a0, a1 = a1, a2 = a2, c = c)
# dat$Y <- (runif(nrow(dat)) < dat$p) + 0
# dat$Group <- (dat$biomarker >= cp) + 0
#
# library(ggplot2)
# p <- ggplot(data = dat, aes(x = biomarker, y = p)) +
#   geom_point(aes(col = as.factor(Group)))
# p
#
# res <- dat %>%
#  group_by(Group) %>%
#  summarise(meanP = mean(p), meanY = mean(Y))
# res

## Bayesian Model
BHM <- function(datList) {
  ## ---------------------------------------------
  ## -- dat
  ## --  (1) K = #. of groups
  ## --  (2) y = #. of responses
  ## --  (3) n = #. of patients
  ## --  (4) epsilon = value of hayperparameter
  ## --------------------------------------------
  Model <- "
  model{
    for (k in 1:K) {
      y[k] ~ dbin(p[k], n[k])
      p[k] ~ dbeta(alpha, beta)
    }
    alpha ~ dgamma(epsilon, epsilon)
    beta ~ dgamma(epsilon, epsilon)
  }
  "

  Fit <- jags(
    model.file = textConnection(Model),
    data = datList,
    parameters.to.save = c("p", "alpha", "beta"),
    n.iter = 20000,
    n.burnin = 10000,
    n.chains = 1, n.thin = 1, progress.bar = "none")
  PoP <- Fit$BUGSoutput$sims.list$p

  return(PoP)
}


simNumeric <- function(nSim = nSim, n,
                       bMin, bMax,
                       a0, a1, a2, c,
                       p0, p1) {
  #nSim <- 1
  set.seed(nSim)
  ## Generate data
  dat <- GenTrialData(n = n, bMin = bMin, bMax = bMax,
                      a0 = a0, a1 = a1, a2 = a2, c = c)
  dat$y <- (runif(nrow(dat)) <= dat$p) + 0

  #postscript(paste(getwd(),"/Result/","DensityPlots_",nSim,".eps", sep = ""),
  #           width = 9, height = 9, onefile = FALSE)
  #par(mfrow = c(3,3))
  ## PoP(p1 - p2|Data)
  cp <- quantile(dat$biomarker, prob = seq(0.1, 0.9, by = 0.1))
  ## cp <- seq(0.5, 4.5, by = 0.5)
  PoP1 <- PoP0 <- PoPDiff <- HPDDiff <- KLDDiff <- NULL
  #PLoss00 <- PLoss01 <- PLoss11 <- NULL
  for (i in 1:length(cp)) {
    dat$Group <- (dat$biomarker >= cp[i]) + 0
    #print(cp[i])
    #dat$Group <- (dat$biomarker >= 2.5) + 0
    Y <- c(sum(dat[dat$Group == 0,]$y),
           sum(dat[dat$Group == 1,]$y))
    N <- c(length(dat[dat$Group == 0,]$y),
           length(dat[dat$Group == 1,]$y))
    cat("# of responses: y(+) = ",Y[2], " and y(-) = ", Y[1], "\n")
    cat("# of patients: n(+) = ",N[2], " and n(-) = ", N[1], "\n")
    cat("PoP(p+|D) = ", round((0.01 + Y[2])/(0.02 + N[2]),4),
        " and PoP(p-|D) = ", round((0.01 + Y[1])/(0.02 + N[1]),4), "\n")
    #print((0.01 + Y[2])/(0.02 + N[2]))
    #print((0.01 + Y[1])/(0.02 + N[1]))

    if (sum(Y == N) <= 0) {
      # datList <- list(K = 2,
      #                 y = c(sum(dat[dat$Group == 0,]$y),
      #                       sum(dat[dat$Group == 1,]$y)),
      #                 n = c(length(dat[dat$Group == 0,]$y),
      #                       length(dat[dat$Group == 1,]$y)),
      #                 epsilon = 0.01)
      # PoP <- BHM(datList = datList)

      y1 <- sum(dat[dat$Group == 1,]$y)
      n1 <- length(dat[dat$Group == 1,]$y)

      y0 <- sum(dat[dat$Group == 0,]$y)
      n0 <- length(dat[dat$Group == 0,]$y)

      PoP <- as.matrix(cbind(rbeta(20000, 0.01 + y0, 0.01 + n0 - y0),
                             rbeta(20000, 0.01 + y1, 0.01 + n1 - y1)))


      x <- seq(0,1, length.out = 1000)
      betaPDF <- function(alpha, beta, x) {
        return( round( x^(alpha - 1)*(1 - x)^(beta - 1)/beta(alpha, beta), 6) )
      }
      pop1 <- lapply(1:length(x),
                    function(i) betaPDF(alpha = 0.01 + y1, beta = 0.01 + n1 - y1, x[i]))
      pop1 <- unlist(pop1)
      pop0 <- lapply(1:length(x),
                     function(i) betaPDF(alpha = 0.01 + y0, beta = 0.01 + n0 - y0, x[i]))
      pop0 <- unlist(pop0)


      #plot(x = x, y = pop0, type = "l",
      #     xlab = "p", ylab = "Density",
      #     col = "blue", ylim = c(0, ceiling(max(c(pop1,pop0)))),
      #     main = paste("Cutoff value at ", names(cp)[i], " percentile", sep = ""))
      #lines(x = x, y = pop1, type = "l", col = "red")
      #legend("topright",
      #       legend = c(expression(p["-"]),expression(p["+"])),
      #       col = c("blue","red"),
      #       lty = 1, bty = "n")


      ## PoP
      PoP1 <- c(PoP1, (0.01 + y1)/(0.02 + n1))
      PoP0 <- c(PoP0, (0.01 + y0)/(0.02 + n0))

      ## PoP Difference
      PoPDiff <- c(PoPDiff, mean(PoP[,2] > PoP[,1]))

      #PLoss00 <- c(PLoss00, sqrt(n1/(n1 + n0))*mean(PoP[,2] >= p0) + sqrt(n0/(n1 + n0))*mean(PoP[,1] < p0))
      #PLoss01 <- c(PLoss01, sqrt(n1/(n1 + n0))*mean(PoP[,2] >= p0) + sqrt(n0/(n1 + n0))*mean(PoP[,1] < p1))
      #PLoss11 <- c(PLoss11, sqrt(n1/(n1 + n0))*mean(PoP[,2] >= p1) + sqrt(n0/(n1 + n0))*mean(PoP[,1] < p1))
      ## HDI
      HPDplus <- emp.hpd(PoP[,2], conf = 0.95)
      HPDminu <- emp.hpd(PoP[,1], conf = 0.95)
      #HPDplus <- hpd(qbeta, shape1 = 0.01 + y1, shape2 = 0.01 + n1 - y1)
      #HPDminu <- hpd(qbeta, shape1 = 0.01 + y0, shape2 = 0.01 + n0 - y0)
      if (HPDplus[1] > HPDminu[2]) {
        Diff <- 0
      } else if (HPDplus[1] <= HPDminu[2] & HPDplus[1] >= HPDminu[1] & HPDminu[2] <= HPDplus[2]){
        Diff <- abs(HPDminu[2] - HPDplus[1])
      } else if (HPDplus[2] <= HPDminu[2] & HPDplus[2] >= HPDminu[1] & HPDplus[1] <= HPDminu[1]) {
        Diff <- abs(HPDplus[2] - HPDminu[1])
      } else if (HPDplus[2] <= HPDminu[2] & HPDplus[1] >= HPDminu[1]) {
        Diff <- abs(HPDplus[2] - HPDplus[1])
      } else if (HPDplus[1] <= HPDminu[1] & HPDplus[2] >= HPDminu[2]) {
        Diff <- abs(HPDminu[2] - HPDminu[1])
      }
      HPDDiff <- c(HPDDiff,  Diff)

      ## HLD difference
      KLDD <- KLx.divergence(PoP[,2],PoP[,1])
      KLDD <- KLDD[which(KLDD != "Inf" & !is.na(KLDD) & KLDD != "-Inf")]
      if (length(KLDD) > 0) {
        KLDDiff <- c(KLDDiff, mean(KLDD))
      } else {
        KLDDiff <- c(KLDDiff, 0)
      }
      #stop(paste("The optimal cutoff value ", format(cp[i],4), " .", sep = ""))
    } else {
      break
    }
  }
  #dev.off()

  if (length(KLDDiff) == length(cp)) {
    ResDiff <- data.frame(
      CutOff = cp,
      PoP1 = PoP1, PoP0 = PoP0,
      #PLoss00 = PLoss00,
      #PLOss01 = PLoss01,
      #PLoss11 = PLoss11,
      PoPDiff = PoPDiff,
      HPDDiff = HPDDiff,
      KLDDiff = KLDDiff
    )
   } else {
     ResDiff <- data.frame(
       CutOff = cp[1:length(KLDDiff)],
       PoP1 = PoP1, PoP0 = PoP0,
       #PLoss00 = PLoss00,
       #PLOss01 = PLoss01,
       #PLoss11 = PLoss11,
       PoPDiff = PoPDiff,
       HPDDiff = HPDDiff,
       KLDDiff = KLDDiff
     )
   }

  #View(ResDiff)
  names(ResDiff) <- c("CutOff",
                      "PoP(P+|D)", "PoP(P-|D)",
                      #"PLoss00","PLoss01","PLoss02",
                      "Pr(P+ > P-|D)",
                      "HPD(P+|D) - HPD(P-|D)",
                      "KLD(P+ - P-|D)")
  #print(ResDiff)
  return(ResDiff)
}

# ResDiff <- simNumeric(nSim = 1, n = n, bMin = bMin, bMax = bMax,
#                       a0 = a0, a1 = a1, a2 = a2, c = c, p0 = p0, p1 = p1)

## Simulation
# nSimulation <- 1000
# Res1 <- mclapply(1:nSimulation,
#                 function(nSim)   simNumeric(nSim = nSim, n = n,
#                                             bMin = bMin, bMax = bMax,
#                                             a0 = a0, a1 = a1, a2 = a2,
#                                             c = c, p0 = p0, p1 = p1),
#                 mc.cores = 4)

nSimulation <- 10
Res1 <- do.call(
  "rbind",
  lapply(1:nSimulation,
         function(nSim)   simNumeric(nSim = nSim, n = n,
                                     bMin = bMin, bMax = bMax,
                                     a0 = a0, a1 = a1, a2 = a2,
                                     c = c, p0 = p0, p1 = p1))
)


ResMat <- matrix(0, nrow = 9, ncol = 7)
Effect <- rep(0,9)
for (i in 1:nSimulation) {
  if (class(Res1[[i]]) != "try-error") {
    res <- as.matrix(Res1[[i]])
    Effect <- Effect + c(rep(1, nrow(res)), rep(0, 9 - nrow(res)))
    ResMat[1:nrow(res),1:ncol(res)] <- ResMat[1:nrow(res),1:ncol(res)] + res
  }
}
ResMat <- as.data.frame(ResMat/Effect)
names(ResMat) <- c("CutOff",
                   "PoP(P+|D)", "PoP(P-|D)",
                   #"PLoss00","PLoss01","PLoss11",
                   "Pr(P+ > P-|D)", "HPD(P+|D) - HPD(P-|D)",
                   "KLD(P+ - P-|D)")
View(round(ResMat,6))
#write.csv(ResMat[,-c(4,5,6)],paste(getwd(),"/Result/","Res_",n,"_NewUpdate.csv", sep = ""))



##
alpha <- 0.01
beta <- 60.01
x <- seq(0,1, length.out = 1000)
betaPDF <- function(alpha, beta, x) {
  return( round( x^(alpha - 1)*(1 - x)^(beta - 1)/beta(alpha, beta), 6) )
}
pop <- lapply(1:length(x),
              function(i) betaPDF(alpha = alpha, beta = beta, x[i]))
pop <- unlist(pop)


plot(x = x, y = pop, type = "l", xlab = "p", ylab = "density")


