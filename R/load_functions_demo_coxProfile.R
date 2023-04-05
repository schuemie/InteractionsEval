library(survival)


##simulate data###
createSetting2DCox <- function (nSites = 5, n = 10000, treatedFraction = 0.2, maleFraction = 0.2, CensorHazard=0.001,
                             minBackgroundHazard = 2e-07, maxBackgroundHazard = 2e-05, 
                             trtFBeta = log(2), trtMBeta = log(2), MBeta = log(2),
                             randomEffectSd = 0, nStrata = 1) {
  expand <- function(x) {
    if (length(x) == 1) {
      return(rep(x, nSites))
    }
    else {
      return(x)
    }
  }
  settings <- list(nSites = nSites, n = expand(n), treatedFraction = expand(treatedFraction), 
                   maleFraction = expand(maleFraction),
                   CensorHazard =CensorHazard,
                   minBackgroundHazard = expand(minBackgroundHazard), 
                   maxBackgroundHazard = expand(maxBackgroundHazard), trtFBeta = expand(trtFBeta),
                   trtMBeta  = expand(trtMBeta ), MBeta = expand(MBeta),
                   randomEffectSd = randomEffectSd, nStrata = expand(nStrata))
  class(settings) <- "simulationSettings"
  return(settings)
}


simulatePopulations2DCox <- function (settings = createSimulationSettings()) 
{
  stopifnot(class(settings) == "simulationSettings")
  #trtBetas <- rnorm(n = settings$nSites, mean = settings$trtBeta, sd = settings$randomEffectSd)
  trtFBetas <- settings$trtFBeta
  #maleBetas <- rnorm(n = settings$nSites, mean = settings$maleBeta, sd = settings$randomEffectSd)
  trtMBetas<-settings$trtMBeta
  #interBetas <- rnorm(n = settings$nSites, mean = settings$interBeta, sd = settings$randomEffectSd)
  
  simulateSite <- function(i) {
    population <- data.frame(rowId = 1:settings$n[i],siteID=i, stratumID = round(runif(settings$n[i], min = 1, max = settings$nStrata[i])), 
                             y = 0, trt = as.numeric(runif(settings$n[i]) < settings$treatedFraction[i]),
                             male = as.numeric(runif(settings$n[i]) < settings$maleFraction[i]))
    strataBackgroundHazard <- runif(settings$nStrata[i], 
                                    min = settings$minBackgroundHazard[i], max = settings$maxBackgroundHazard[i])
    MBetas <- rnorm(n = settings$nStrata[i], mean = settings$MBeta, sd = settings$randomEffectSd)
    population$hazard <- strataBackgroundHazard[population$stratumID]
    #oldTotalHazard <- sum(population$hazard)
    population$inter<-population$trt*population$male
    #covariates <- cbind(population$trt,population$male, population$trt*population$male)
    hazardratios <- exp(population$trt*(1-population$male)*trtFBetas[i]+population$inter*trtMBetas[i]
                        +population$male*MBetas[population$stratumID])
    population$hazard <- population$hazard * hazardratios
    
    
    #newTotalHazard <- sum(population$hazard)
    #population$hazard <- population$hazard * oldTotalHazard/newTotalHazard
    population$timeToOutcome <- 1 + round(rexp(n = settings$n[i], 
                                               population$hazard))
    #CensorHazard<-2.9*quantile(population$hazard,settings$CensorFraction)
    population$timeToCensor <- 1 + round(rexp(n = settings$n[i], 
                                              settings$CensorHazard))
    population$time <- population$timeToOutcome
    population$x1<-population$trt*(1-population$male)
    population$x2<-population$inter
    population$z1<-population$male
    population$Treatment<-population$trt
    idx <- population$timeToCensor < population$timeToOutcome
    population$time[idx] <- population$timeToCensor[idx]
    population$status <- as.integer(!idx)
    return(population[, c("time","status","x1","x2","z1","siteID", "stratumID", "Treatment","male")])
  }
  simulation <- lapply(1:settings$nSites,simulateSite)
  attr(simulation, "simulationSettings") <- settings
  attr(simulation, "trtFBetas") <- trtFBetas
  attr(simulation, "trtMBetas") <- trtMBetas
  class(simulation) <- "simulation"
  return(simulation)
}





# Pooled analysis
nlogLp_poolCox <- function(eta,ppdata){
  Lp<-rep(0)
  SiteID<-ppdata$siteID
  J<- max(SiteID)
  for (j in 1:J){
    spdata<-ppdata[SiteID==j,]
    StratumId<-spdata$stratumID
    K<-max(StratumId)
    for(k in 1:K){
      ipdata <- spdata[StratumId==k,]
      ipdata<-ipdata[order(-ipdata$time),]
      y_time<-ipdata$time
      y_status<-ipdata$status
      index_1<-which(y_status==TRUE)
      nlogLp <- function(eta, ipdata){
        x1<-ipdata$x1
        x2<-ipdata$x2
        x3<-ipdata$z1
        X <- cbind(x1,x2,x3)
        fit_j<-coxph.fit(x=cbind(X[,3]), y=Surv(y_time, y_status),strata = NULL,offset = X[,1:2]%*%eta,
                         control =coxph.control(),rownames = NULL,method="breslow")
        gamma_fit<-fit_j$coefficients
        beta<-c(eta,gamma_fit)
        y_mu<-X%*%beta
        expy<-exp(y_mu)
        
        if(length(index_1)>0){
          s_expy<-cumsum(expy)[index_1]
          nlogLp <- -sum(y_mu[index_1])+sum(log(s_expy))
        }else{
          nlogLp<-0
        }
        return(nlogLp)
      }
      Lp <- nlogLp(eta,ipdata)+Lp
      #print(nlogLp(eta))
    }
  }
  return(Lp)
}


estimatePoolCox<- function(ppdata){
  eta.pool <- optim(c(1,1),nlogLp_poolCox,ppdata = ppdata)$par
  return(list(est = c(eta.pool)))
}





####meta###
fit.barCox <- function(ppdata){  ## use glm
  J<- max(ppdata$siteID)
  ehat <- matrix(NA,nrow = J,ncol=2)
  Vhat <-array(NA,c(J,2,2))
  for (j in 1:J){
    ipdata <- ppdata[ppdata$siteID==j,]
    x1<-ipdata$x1
    x2<-ipdata$x2
    x3<-ipdata$z1
    StratumId <- ipdata$stratumID
    new.stratumID <- sapply(1:max(StratumId),function(i) ifelse(StratumId==i,1,0))
    X <- cbind(x1,x2,x3*new.stratumID)
    #X <- cbind(x1,x2,x3)
    y_time<-ipdata$time
    y_status<-ipdata$status
    fit_j<-coxph.fit(x=X, y=Surv(y_time, y_status),strata = NULL,
                     control =coxph.control(),rownames = NULL, method="breslow")
    ehat[j,] = fit_j$coefficients[1:2]
    Vhat[j,,] = fit_j$var[1:2,1:2]
  }
  return(list(ehat=ehat,Vhat=Vhat))
}

estimateMetaCox <- function(ppdata){
  fit.local<-fit.barCox(ppdata)
  # remove studies with null
  #ValidSite<-which(rowSums(is.na(fit.local$ehat))==0)
  ValidSite<-which(fit.local$Vhat[,1,1]>0  & fit.local$Vhat[,2,2]>0 )
  eta.local<- fit.local$ehat[ValidSite,]
  vcov.local <- fit.local$Vhat[ValidSite,,]
  eMeta <- rep(0)
  wtMeta <- matrix(0,2,2)
  for(j in 1:length(ValidSite)){
    eMeta<- eMeta+solve(vcov.local[j,,])%*%eta.local[j,]
    wtMeta<- wtMeta+solve(vcov.local[j,,])
  }
  eMeta <- solve(wtMeta)%*%eMeta
  vMeta<- solve(wtMeta)
  return(list(est = c(eMeta),cov = vMeta))
}



###pade###

logLp_derivCox <- function(ipdata,eta){
  ipdata<-ipdata[order(-ipdata$time),]
  y_time<-ipdata$time
  y_status<-ipdata$status
  index_1<-which(y_status==TRUE)
  x1<-ipdata$x1
  x2<-ipdata$x2
  x3<-ipdata$z1
  X <- cbind(x1,x2,x3)
  if(length(index_1)>0){
  fit_j<-coxph.fit(x=cbind(X[,3]), y=Surv(y_time, y_status),strata = NULL,offset = X[,1:2]%*%eta,
                   control =coxph.control(),rownames = NULL,method="breslow")
  beta<-c(eta,fit_j$coefficients)
  expy <- exp(X%*%beta)

    s_expy<-cumsum(expy)[index_1]
    expy_z<-cumsum(expy*x3)[index_1]
    expy_x1<-cumsum(expy*x1)[index_1]
    expy_x2<-cumsum(expy*x2)[index_1]
    expy_zz<-cumsum(expy*x3*x3)[index_1]
    expy_zx1<-cumsum(expy*x3*x1 )[index_1]
    expy_zx2<-cumsum(expy*x3*x2 )[index_1]
    
    invl01<-1/(sum( expy_zz/s_expy-expy_z*expy_z/(s_expy^2))) 
    D1_e1 <- -invl01*sum(expy_zx1/s_expy- expy_z*expy_x1/(s_expy^2))
    D1_e2 <- -invl01*sum(expy_zx2/s_expy- expy_z*expy_x2/(s_expy^2))
    
    
    expy_p1<-cumsum(expy*(x1+x3*D1_e1))[index_1]
    expy_p2<-cumsum(expy*(x2+x3*D1_e2))[index_1]
    expy_zp1<-cumsum(expy*x3*(x1+x3*D1_e1))[index_1]
    expy_zp2<-cumsum(expy*x3*(x2+x3*D1_e2))[index_1]
    expy_p1p2<-cumsum(expy*(x1+x3*D1_e1)*(x2+x3*D1_e2))[index_1]
    expy_p1p1<-cumsum(expy*(x1+x3*D1_e1)*(x1+x3*D1_e1))[index_1]
    expy_p2p2<-cumsum(expy*(x2+x3*D1_e2)*(x2+x3*D1_e2))[index_1]
    expy_zp1p2<-cumsum(expy*x3*(x1+x3*D1_e1)*(x2+x3*D1_e2))[index_1]
    expy_zp1p1<-cumsum(expy*x3*(x1+x3*D1_e1)*(x1+x3*D1_e1))[index_1]
    expy_zp2p2<-cumsum(expy*x3*(x2+x3*D1_e2)*(x2+x3*D1_e2))[index_1]
    
    
    D2_e11 <- -invl01*sum(expy_zp1p1/s_expy-(2*expy_zp1*expy_p1+expy_z*expy_p1p1)/(s_expy^2)
                          +2*expy_z*expy_p1*expy_p1/(s_expy^3))
    D2_e12 <- -invl01*sum(expy_zp1p2/s_expy-(expy_zp1*expy_p2+expy_zp2*expy_p1+expy_z*expy_p1p2)/(s_expy^2)
                          +2*expy_z*expy_p1*expy_p2/(s_expy^3))
    D2_e22 <- -invl01*sum(expy_zp2p2/s_expy-(2*expy_zp2*expy_p2+expy_z*expy_p2p2)/(s_expy^2)
                          +2*expy_z*expy_p2*expy_p2/(s_expy^3))
    
    
    expy_p11<-cumsum(expy*(x3*D2_e11))[index_1]
    expy_p12<-cumsum(expy*(x3*D2_e12))[index_1]
    expy_p22<-cumsum(expy*(x3*D2_e22))[index_1]
    expy_zp11<-cumsum(expy*x3*(x3*D2_e11))[index_1]
    expy_zp12<-cumsum(expy*x3*(x3*D2_e12))[index_1]
    expy_zp22<-cumsum(expy*x3*(x3*D2_e22))[index_1]
    
    expy_p1p11<-cumsum(expy*(x1+x3*D1_e1)*(x3*D2_e11))[index_1]
    expy_p1p12<-cumsum(expy*(x1+x3*D1_e1)*(x3*D2_e12))[index_1]
    expy_p1p22<-cumsum(expy*(x1+x3*D1_e1)*(x3*D2_e22))[index_1]
    expy_p2p11<-cumsum(expy*(x2+x3*D1_e2)*(x3*D2_e11))[index_1]
    expy_p2p12<-cumsum(expy*(x2+x3*D1_e2)*(x3*D2_e12))[index_1]
    expy_p2p22<-cumsum(expy*(x2+x3*D1_e2)*(x3*D2_e22))[index_1]
    
    expy_zp1p11<-cumsum(expy*x3*(x1+x3*D1_e1)*(x3*D2_e11))[index_1]
    expy_zp1p12<-cumsum(expy*x3*(x1+x3*D1_e1)*(x3*D2_e12))[index_1]
    expy_zp1p22<-cumsum(expy*x3*(x1+x3*D1_e1)*(x3*D2_e22))[index_1]
    expy_zp2p11<-cumsum(expy*x3*(x2+x3*D1_e2)*(x3*D2_e11))[index_1]
    expy_zp2p12<-cumsum(expy*x3*(x2+x3*D1_e2)*(x3*D2_e12))[index_1]
    expy_zp2p22<-cumsum(expy*x3*(x2+x3*D1_e2)*(x3*D2_e22))[index_1]
    
    expy_p1p1p1<-cumsum(expy*(x1+x3*D1_e1)^3)[index_1]
    expy_p1p1p2<-cumsum(expy*(x1+x3*D1_e1)^2*(x2+x3*D1_e2))[index_1]
    expy_p1p2p2<-cumsum(expy*(x1+x3*D1_e1)*(x2+x3*D1_e2)^2)[index_1]
    expy_p2p2p2<-cumsum(expy*(x2+x3*D1_e2)^3)[index_1]
    
    expy_zp1p1p1<-cumsum(expy*x3*(x1+x3*D1_e1)^3)[index_1]
    expy_zp1p1p2<-cumsum(expy*x3*(x1+x3*D1_e1)^2*(x2+x3*D1_e2))[index_1]
    expy_zp1p2p2<-cumsum(expy*x3*(x1+x3*D1_e1)*(x2+x3*D1_e2)^2)[index_1]
    expy_zp2p2p2<-cumsum(expy*x3*(x2+x3*D1_e2)^3)[index_1]
    
    
    D3_e111 <- -invl01*sum((expy_zp1p1p1+3*expy_zp1p11)/s_expy-(3*expy_zp1p1*expy_p1+3*expy_zp1*expy_p1p1
                                                                +3*expy_zp11*expy_p1+3*expy_p11*expy_zp1+3*expy_z*expy_p1p11+expy_z*expy_p1p1p1)/(s_expy^2)+
                             (6*expy_zp1*expy_p1^2+6*expy_z*expy_p1p1*expy_p1+6*expy_z*expy_p1*expy_p11)/(s_expy^3)-
                             6*(expy_z*expy_p1^3)/(s_expy^4))
    
    
    D3_e112 <- -invl01*sum((expy_zp1p1p2+2*expy_zp1p12+expy_zp2p11)/s_expy-(2*expy_zp1p2*expy_p1+expy_zp1p1*expy_p2
                                                                            +2*expy_zp1*expy_p1p2+expy_zp2*expy_p1p1+2*expy_zp12*expy_p1+expy_zp11*expy_p2+2*expy_p12*expy_zp1+
                                                                              expy_p11*expy_zp2+2*expy_z*expy_p1p12+expy_z*expy_p2p11+expy_z*expy_p1p1p2)/(s_expy^2)+
                             2*(2*expy_zp1*expy_p1*expy_p2+expy_zp2*expy_p1*expy_p1+2*expy_z*expy_p1p2*expy_p1+expy_z*expy_p1p1*expy_p2
                                +2*expy_z*expy_p1*expy_p12+expy_z*expy_p2*expy_p11)/(s_expy^3)-
                             6*(expy_z*expy_p1*expy_p1*expy_p2)/(s_expy^4))
    
    D3_e122 <- -invl01*sum((expy_zp1p2p2+2*expy_zp2p12+expy_zp1p22)/s_expy-(2*expy_zp1p2*expy_p2+expy_zp2p2*expy_p1
                                                                            +2*expy_zp2*expy_p1p2+expy_zp1*expy_p2p2+2*expy_zp12*expy_p2+expy_zp22*expy_p1+2*expy_p12*expy_zp2+
                                                                              expy_p22*expy_zp1+2*expy_z*expy_p2p12+expy_z*expy_p1p22+expy_z*expy_p1p2p2)/(s_expy^2)+
                             2*(2*expy_zp2*expy_p1*expy_p2+expy_zp1*expy_p2*expy_p2+2*expy_z*expy_p1p2*expy_p2+expy_z*expy_p2p2*expy_p1
                                +2*expy_z*expy_p2*expy_p12+expy_z*expy_p1*expy_p22)/(s_expy^3)-
                             6*(expy_z*expy_p1*expy_p2*expy_p2)/(s_expy^4))
    
    D3_e222 <- -invl01*sum((expy_zp2p2p2+3*expy_zp2p22)/s_expy-(3*expy_zp2p2*expy_p2+3*expy_zp2*expy_p2p2
                                                                +3*expy_zp22*expy_p2+3*expy_p22*expy_zp2+3*expy_z*expy_p2p22+expy_z*expy_p2p2p2)/(s_expy^2)+
                             (6*expy_zp2*expy_p2^2+6*expy_z*expy_p2p2*expy_p2+6*expy_z*expy_p2*expy_p22)/(s_expy^3)-
                             6*(expy_z*expy_p2^3)/(s_expy^4))
    
    logL <- sum((X%*%beta)[index_1])-sum(log(s_expy))
    logL_D1 <-c(sum(x1[index_1]+x3[index_1]*D1_e1-expy_p1/s_expy),sum(x2[index_1]+x3[index_1]*D1_e2-expy_p2/s_expy))
    logL_D2_11 <-sum(x3[index_1]*D2_e11-(expy_p1p1+expy_p11)/s_expy+expy_p1*expy_p1/(s_expy^2)) 
    logL_D2_12 <-sum(x3[index_1]*D2_e12-(expy_p1p2+expy_p12)/s_expy+expy_p1*expy_p2/(s_expy^2)) 
    logL_D2_22 <-sum(x3[index_1]*D2_e22-(expy_p2p2+expy_p22)/s_expy+expy_p2*expy_p2/(s_expy^2))
    logL_D2 <- matrix(c(logL_D2_11,logL_D2_12,logL_D2_12,logL_D2_22),nrow=2,ncol=2)
    
    logL_D3_111 <-sum(x3[index_1]*D3_e111-expy_z*D3_e111/s_expy -(expy_p1p1p1+3*expy_p1p11)/(s_expy)+(3*expy_p1*expy_p1p1+
                                                                                                        3*expy_p1*expy_p11)/(s_expy^2) -2*(expy_p1^3)/(s_expy^3) )
    logL_D3_112 <- sum(x3[index_1]*D3_e112-expy_z*D3_e112/s_expy-(expy_p1p1p2+2*expy_p1p12+expy_p2p11)/(s_expy)+(2*expy_p1*expy_p1p2+
                                                                                                                   expy_p2*expy_p1p1+2*expy_p1*expy_p12+expy_p2*expy_p11)/(s_expy^2) -2*(expy_p1^2*expy_p2)/(s_expy^3) )
    logL_D3_122 <- sum(x3[index_1]*D3_e122-expy_z*D3_e122/s_expy -(expy_p1p2p2+2*expy_p2p12+expy_p1p22)/(s_expy)+(2*expy_p2*expy_p1p2+
                                                                                                                    expy_p1*expy_p2p2+2*expy_p2*expy_p12+expy_p1*expy_p22)/(s_expy^2) -2*(expy_p2^2*expy_p1)/(s_expy^3))
    logL_D3_222 <-sum(x3[index_1]*D3_e222-expy_z*D3_e222/s_expy -(expy_p2p2p2+3*expy_p2p22)/(s_expy)+(3*expy_p2*expy_p2p2+
                                                                                                        3*expy_p2*expy_p22)/(s_expy^2) -2*(expy_p2^3)/(s_expy^3) )
  }else{
    logL<-0
    logL_D1<-c(0,0)
    logL_D2<-c(0,0)
    logL_D2<-matrix(0,2,2)
    logL_D3_111<-0
    logL_D3_222<-0
    logL_D3_112<-0
    logL_D3_122<-0
  }
  return(list(logLp=logL,logLp_D1=logL_D1,logLp_D2=logL_D2,logLp_D3=c(logL_D3_111,logL_D3_222,logL_D3_112,logL_D3_122)))
}


GetLocalDerivCox <- function(ebar,ipdata){
  StratumId<-ipdata$stratumID
  K<- max(StratumId)
  logLp <- 0
  logLp_D1 <- rep(0,2)
  logLp_D2 <- matrix(0,2,2)
  logLp_D3 <- rep(0,4)
  for (k in 1:K){
    ipdata.strata <- ipdata[StratumId==k,]
    logLp_derivCox_j<-logLp_derivCox(ipdata.strata,ebar)
    logLp <- logLp_derivCox_j$logLp+logLp
    logLp_D1 <- logLp_derivCox_j$logLp_D1+logLp_D1
    logLp_D2 <- logLp_derivCox_j$logLp_D2+logLp_D2
    logLp_D3 <- logLp_derivCox_j$logLp_D3+logLp_D3
  }
  return(list(logLp=logLp,logLp_D1=logLp_D1,logLp_D2=logLp_D2,logLp_D3=logLp_D3))
}




GetGlobalDerivCox <- function(ebar, ppdata){
  SiteID<-ppdata$siteID
  J<- max(SiteID)
  logLp <- 0
  logLp_D1 <- rep(0,2)
  logLp_D2 <- matrix(0,2,2)
  logLp_D3 <- rep(0,4)
  for (j in 1:J){
    spdata<-ppdata[SiteID==j,]
    StratumId<-spdata$stratumID
    K<-max(StratumId)
    for(k in 1:K){
      ipdata <- spdata[StratumId==k,]
      logLp_derivCox_j<-logLp_derivCox(ipdata,ebar)
      logLp <- logLp_derivCox_j$logLp+logLp
      logLp_D1 <- logLp_derivCox_j$logLp_D1+logLp_D1
      logLp_D2 <- logLp_derivCox_j$logLp_D2+logLp_D2
      logLp_D3 <- logLp_derivCox_j$logLp_D3+logLp_D3
    }
  }
  return(list(logLp=logLp,logLp_D1=logLp_D1,logLp_D2=logLp_D2,logLp_D3=logLp_D3))
}

EstPadeCoefCox<- function(deriv,Midx = rbind(c(0,0),c(1,0),c(0,1),c(1,1),c(2,0),c(0,2)),
                       Nidx = rbind(c(0,0),c(1,0),c(0,1),c(1,1),c(2,0),c(0,2)),
                       Eidx = rbind(Midx,c(3,0),c(0,3))){
  c0 <- deriv$logLp
  c1 <- deriv$logLp_D1
  c2 <- deriv$logLp_D2/2
  Esub_idx<-Eidx[c((dim(Midx)[1]+1):(dim(Eidx)[1])),]
  Esub_N<-array(0,c(dim(Esub_idx),dim(Nidx)[1]-1))
  
  A <- matrix(0,dim(Esub_idx)[1],(dim(Nidx)[1]-1))
  for(i in 1:dim(Esub_idx)[1]){
    for(j in 1:(dim(Nidx)[1]-1)){
      EN_idx <- Esub_idx[i,]-Nidx[j+1,]
      if(min(EN_idx)<0) A[i,j]=0
      else if (sum(EN_idx)==0) A[i,j]<-c0
      else if (sum(EN_idx)==1) A[i,j]<-c1[which(EN_idx==1)]
      else if (any(EN_idx==2)) A[i,j]<-c2[which(EN_idx==2),which(EN_idx==2)]
      else{
        idx2 <- which(EN_idx==1)
        A[i,j]<-2*c2[idx2[1],idx2[2]]
      }
    }
  }
  
  A2 = matrix(0,dim(Midx)[1],dim(Nidx)[1])
  for(i in 1:dim(Midx)[1]){
    for(j in 1:(dim(Nidx)[1])){
      MN_idx <- Midx[i,]-Nidx[j,]
      if(min(MN_idx)<0) A2[i,j]=0
      else if (sum(MN_idx)==0) A2[i,j]=c0
      else if (sum(MN_idx)==1) A2[i,j]=c1[which(MN_idx==1)]
      else if (any(MN_idx==2)) A2[i,j]=c2[which(MN_idx==2),which(MN_idx==2)]
      else{
        idx2 <- which(MN_idx==1)
        A2[i,j]=2*c2[idx2[1],idx2[2]]
      }
    }
  }
  A00 <- c(deriv$logLp_D3[1:2]/6)
  A <- cbind(A00,A)
  
  A_decomp <- eigen(t(A)%*%A)
  min_idx <- which(abs(A_decomp$vectors[1,])>1e-4 & abs(A_decomp$values)<1e-8)
  min_idx <- min_idx[which.min(abs(A_decomp$values[min_idx]))]
  vv<- A_decomp$vectors[,min_idx]
  
  vvQ<-vv
  vvP <- A2%*%vvQ
  return(list(Pcoef = c(vvP),Qcoef = c(vvQ)))
}


PadeEstCRCox<- function(ebar,PadeCoef,alpha=0.05){
  num_vec<-c(PadeCoef$Pcoef)
  denom_vec<-c(PadeCoef$Qcoef)
  Pade<-function(eta){ # pade function to approximate negative log-likelihood
    x <- eta-ebar
    terms <- c(1,x,x[1]*x[2],x^2)
    Pm<- sum(terms*num_vec)
    Qn <- sum(terms*denom_vec)
    pade <- -Pm/Qn
    return(pade)
  }
  eta.pade<-optim(ebar,Pade,method = "L-BFGS-B",
                  lower = ebar-log(10), upper = ebar+log(10))
  PadeEst <- eta.pade$par
  
  # CR
  coef.CR<-num_vec-denom_vec*c(-eta.pade$value-qchisq(1-alpha,2)/2)
  
  return(list(Est=PadeEst,CR=coef.CR))
}


## function to plot CR

ellipse.plot.metaCox <- function(mua, mub, sa, sb, rho, n = NULL, alpha=0.05){
  pi <- 3.14159265359
  tt<-seq(0, 2*pi, length=1000)
  #CC<-sqrt(qchisq(1-alpha, 2))
  #CC<-sqrt(2*qf(1-alpha, 2, n-2))
  CC <- qnorm(1-alpha/2)
  ua <<- mua + sa*CC*cos(tt)
  ub <<- mub + sb*CC*cos(tt+acos(rho))
  return(data.frame(x = ua,y = ub))
}

# ellipse.plot.padeCox0Cox: old code to draw ellipse
# source: https://stackoverflow.com/questions/41820683/how-to-plot-ellipse-given-a-general-equation-in-r
ellipse.plot.padeCox0Cox <- function(coef, ebar){
  
  A = coef[5]
  C = coef[6]
  B = coef[4]
  D = coef[2]
  E = coef[3]
  FF =coef[1]
  pi <- 3.14159265359
  M0 <- matrix(c(FF, D/2, E/2, D/2, A, B/2, E/2, B/2, C), 3L,
               3L)
  M <- matrix(c(A, B/2, B/2, C), 2L, 2L)
  lambda <- eigen(M, symmetric = TRUE)$values
  detM0 <- det(M0)
  detM <- det(M)
  
  a <- sqrt(-det(M0)/(det(M)*lambda[1]))
  b <- sqrt(-det(M0)/(det(M)*lambda[2]))
  
  
  x <- B * E - 2 * C * D
  y <- B * D - 2 * A * E
  
  phi <- if (is.nan(B/(A - C))) {
    0
  }else {
    if (abs(C) > abs(A))
      atan(B/(A - C))/2
    else (pi/2 - atan(-B/(A - C))/2)
  }
  xc <- x/(4 * A * C - B * B)
  yc <- y/(4 * A * C - B * B)
  
  
  t <- seq(0, 2*pi, length = 1000)
  x <- xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi)
  y <- yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi)
  return(data.frame(x = x+ebar[1],y = y+ebar[2]))
}


ellipse.plot.padeCox <- function(coef, ebar){
  
  a = coef[5]
  b = coef[4]
  c = coef[6]
  d = coef[2]
  e = coef[3]
  f = coef[1]
  
  pi <- 3.14159265359
  t <- seq(0, 2*pi, length = 1000)
  
  y.dis <- c-b^2/4/a
  yc <- (e-d*b/2/a)/2/y.dis
  yc2 <- yc^2*y.dis-f+d^2/4/a
  
  y <- sqrt(yc2/y.dis)*sin(t)-yc
  x <- sqrt(yc2/a)*cos(t)-d/2/a-b/2/a*y
  return(data.frame(x = x+ebar[1],y = y+ebar[2]))
}
