
SimulateData<- function(N,nSites,nStrata,MaleFraction,TreatedFraction,eta, prevalence = 0.2, followup = 2){
  
  pi <- runif(nSites*nStrata, 0.1,0.9)
  pi <- pi/sum(pi)
  # set nuisance parameters based on prevalence
  gamma2 <- c(0,0.5) # mean and sd of second nuisance parameter
  
  SitestratumIDx <- rmultinom(N,1,pi)
  SitestratumID <- max.col(t(SitestratumIDx),"first")
  siteID <- ceiling(SitestratumID/nStrata)
  stratumID <- SitestratumID-nStrata*(siteID-1)
  
  Male <- rbinom(N,1,MaleFraction)
  Treatment <- rbinom(N,1,TreatedFraction)
  period <- ceiling(rexp(N,1/followup))
  #period <- rep(2,N)
  X <- cbind(Treatment*(1-Male),Treatment*Male,rep(1,N),Male) # first two are parameter of interest (denoted as eta); last two gamma
  colnames(X) <- c("x1", "x2","z1","z2")
  
  J <- length(pi)
  gamma <-matrix(NA,nrow = J,ncol=2)
  theta <-matrix(NA,nrow = J,ncol=4)
  
  y <- rep(NA,N)
  for (j in 1:J){
    nj <- sum(siteID==j)
    gamma[j,2] <- rnorm(1,gamma2[1],gamma2[2])
    gamma[j,1] <- log(prevalence/followup)-eta[1]*TreatedFraction*(1-MaleFraction)-eta[2]*TreatedFraction-gamma[j,2]*MaleFraction
    theta[j,] <- c(eta,gamma[j,])
    mu_j <- exp(X[siteID==j,]%*%theta[j,])
    period_j <- period[siteID==j]
    yj <- rpois(nj, lambda=mu_j*period_j)
    y[siteID==j]=yj
  }
  ppdata <- data.frame(y,X,period,siteID,stratumID, Treatment, Male)
  return(list(ppdata=ppdata, theta = theta))
}

# Pooled analysis
nlogLp_pool <- function(eta,ppdata){
  Lp<-rep(0)
  nStrata <- max(ppdata$stratumID)
  SitestratumID <- (ppdata$siteID-1)*nStrata+ppdata$stratumID
  J<- max(SitestratumID)
  for (j in 1:J){
    ipdata <- ppdata[SitestratumID==j,]
    nlogLp <- function(eta){
      eta0<-eta
      X <- as.matrix(ipdata%>%select("x1","x2","z1","z2"))
      y <- ipdata$y
      period <- ipdata$period
      fit_j<-glm(y ~ -1+X[,3]+X[,4], offset = X[,1:2]%*%eta0+log(period),  family = poisson(link = log))
      gamma_fit<-fit_j$coefficients
      nlogLp <- -sum(y*(X%*%c(eta0,gamma_fit)+log(period))-period*c(exp(X%*%c(eta0,gamma_fit))))
      return(nlogLp)
    }
    Lp <- nlogLp(eta)+Lp
    #print(nlogLp(eta))
  }
  return(Lp)
}
estimatePool<- function(ppdata){
  eta.pool <- optim(c(1,1),nlogLp_pool,ppdata = ppdata)$par
  return(list(est = c(eta.pool)))
}


# Meta analysis
fit.bar <- function(ppdata){  ## use glm
  J<- max(ppdata$siteID)
  ehat <- matrix(NA,nrow = J,ncol=2)
  Vhat <-array(NA,c(J,2,2))
  for (j in 1:J){
    ipdata <- ppdata[ppdata$siteID==j,]
    X <- as.matrix(ipdata%>%select("x1","x2","z1","z2"))
    y <- ipdata$y
    period <- ipdata$period
    stratumID <- ipdata$stratumID
    new.stratumID <- sapply(1:max(stratumID),function(i) ifelse(stratumID==i,1,0))
    Xnew.strata <- cbind(X[,1:2],X[,3]*new.stratumID,X[,4]*new.stratumID)
    fit_j<-glm(y ~ -1+Xnew.strata, offset = log(period), family = poisson(link = log))
    ehat[j,] = fit_j$coefficients[1:2]
    Vhat[j,,] = vcov(fit_j)[1:2,1:2]
  }
  return(list(ehat=ehat,Vhat=Vhat))
}


estimateMeta <- function(ppdata){
  fit.local<-fit.bar(ppdata)
  # remove studies with null
  ValidSite<-which(rowSums(is.na(fit.local$ehat))==0)
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

## Pade

logLp_deriv <- function(ipdata,theta){
  y <- ipdata$y
  x1 <-ipdata$x1; x2<-ipdata$x2;x3<-ipdata$z1; x4<-ipdata$z2
  period <- ipdata$period
  X <- cbind(x1,x2,x3,x4)
  eXt <- c(exp(X%*%theta))*period
  L01 <- colSums((y-eXt)*cbind(x3,x4))
  L02 <- -colSums(eXt*cbind(x3*x3,x3*x4,x4*x3,x4*x4))
  invL02 <- solve(matrix(L02,nrow=2,ncol=2))
  L11_e1 <- -colSums(eXt*cbind(x3*x1,x4*x1))
  L11_e2 <- -colSums(eXt*cbind(x3*x2,x4*x2))
  D1_e1 <- -invL02%*%L11_e1
  D1_e2 <- -invL02%*%L11_e2
  
  
  W1 <- cbind(x3,x4)%*%D1_e1
  W2 <- cbind(x3,x4)%*%D1_e2
  b11 <- -colSums(c(eXt*(x1+W1)^2)*cbind(x3,x4))
  D2_e11 <- -invL02%*%b11
  b12 <- -colSums(c(eXt*(x1+W1)*(x2+W2))*cbind(x3,x4))
  D2_e12 <- -invL02%*%b12
  D2_e21 <- D2_e12
  b22 <- -colSums(c(eXt*(x2+W2)^2)*cbind(x3,x4))
  D2_e22 <- -invL02%*%b22
  
  W11 <- cbind(x3,x4)%*%D2_e11
  W12 <- cbind(x3,x4)%*%D2_e12
  W21 <- cbind(x3,x4)%*%D2_e21
  W22 <- cbind(x3,x4)%*%D2_e22
  b111 <- -colSums(c(eXt*((x1+W1)^3+3*(x1+W1)*W11))*cbind(x3,x4))
  D3_e111 <- -invL02%*%b111
  b112 <- -colSums(c(eXt*((x2+W2)*(x1+W1)^2+(x2+W2)*W11+2*(x1+W1)*W12))*cbind(x3,x4))
  D3_e112 <- -invL02%*%b112
  b122 <- -colSums(c(eXt*((x2+W2)^2*(x1+W1)+2*(x2+W2)*W12+(x1+W1)*W22))*cbind(x3,x4))
  D3_e122 <- -invL02%*%b122
  b222 <- -colSums(c(eXt*((x2+W2)^3+3*(x2+W2)*W22))*cbind(x3,x4))
  D3_e222 <- -invL02%*%b222
  
  logL <- sum(y*(X%*%theta+log(period))-eXt)
  logL_D1 <- c(sum((y-eXt)*(x1+W1)),sum((y-eXt)*(x2+W2)))
  logL_D2_11 <- sum(-eXt*(x1+W1)^2+(y-eXt)*W11)
  logL_D2_12 <- sum(-eXt*(x2+W2)*(x1+W1)+(y-eXt)*W12)
  logL_D2_22 <- sum(-eXt*(x2+W2)^2+(y-eXt)*W22)
  logL_D2 <- matrix(c(logL_D2_11,logL_D2_12,logL_D2_12,logL_D2_22),nrow=2,ncol=2)
  
  W111 <- cbind(x3,x4)%*%D3_e111
  W112 <- cbind(x3,x4)%*%D3_e112
  W122 <- cbind(x3,x4)%*%D3_e122
  W222 <- cbind(x3,x4)%*%D3_e222
  logL_D3_111 <- sum(-eXt*((x1+W1)^3+3*(x1+W1)*W11)+(y-eXt)*W111)
  logL_D3_112 <- sum(-eXt*((x1+W1)^2*(x2+W2)+2*(x1+W1)*W12+(x2+W2)*W11)+(y-eXt)*W112)
  logL_D3_122 <- sum(-eXt*((x1+W1)*(x2+W2)^2+2*(x2+W2)*W12+(x1+W1)*W22)+(y-eXt)*W122)
  logL_D3_222 <- sum(-eXt*((x2+W2)^3+3*(x2+W2)*W22)+(y-eXt)*W222)
  
  
  return(list(logLp=logL,logLp_D1=logL_D1,logLp_D2=logL_D2,logLp_D3=c(logL_D3_111,logL_D3_222,logL_D3_112,logL_D3_122)))
}

GetLocalDeriv <- function(ebar,ipdata){
  
  J<- max(ipdata$stratumID)
  logLp <- 0
  logLp_D1 <- rep(0,2)
  logLp_D2 <- matrix(0,2,2)
  logLp_D3 <- rep(0,4)
  for (j in 1:J){
    ipdata.strata <- ipdata[ipdata$stratumID==j,]
    X <- as.matrix(ipdata.strata%>%select("x1","x2","z1","z2"))
    y <- ipdata.strata$y
    period <- ipdata.strata$period
    
    fit_j<-glm(y ~ -1+X[,3]+X[,4], offset = X[,1:2]%*%ebar+log(period),  family = poisson(link = log))
    gfit_ebar<-fit_j$coefficients
    thetaj <- c(ebar,gfit_ebar)
    logLp_deriv_j<-logLp_deriv(ipdata.strata,thetaj)
    logLp <- logLp_deriv_j$logLp+logLp
    logLp_D1 <- logLp_deriv_j$logLp_D1+logLp_D1
    logLp_D2 <- logLp_deriv_j$logLp_D2+logLp_D2
    logLp_D3 <- logLp_deriv_j$logLp_D3+logLp_D3
  }
  return(list(logLp=logLp,logLp_D1=logLp_D1,logLp_D2=logLp_D2,logLp_D3=logLp_D3))
}

GetGlobalDeriv <- function(ebar,ppdata){
  J<- max(ppdata$siteID)
  logLp <- 0
  logLp_D1 <- rep(0,2)
  logLp_D2 <- matrix(0,2,2)
  logLp_D3 <- rep(0,4)
  for (j in 1:J){
    ipdata <- ppdata[ppdata$siteID==j,]
    logLp_deriv_j<-GetLocalDeriv(ebar,ipdata)
    logLp <- logLp_deriv_j$logLp+logLp
    logLp_D1 <- logLp_deriv_j$logLp_D1+logLp_D1
    logLp_D2 <- logLp_deriv_j$logLp_D2+logLp_D2
    logLp_D3 <- logLp_deriv_j$logLp_D3+logLp_D3
  }
  return(list(logLp=logLp,logLp_D1=logLp_D1,logLp_D2=logLp_D2,logLp_D3=logLp_D3))
}

EstPadeCoef<- function(deriv,Midx = rbind(c(0,0),c(1,0),c(0,1),c(1,1),c(2,0),c(0,2)),
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

PadeEstCR<- function(ebar,PadeCoef,alpha=0.05){
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

ellipse.plot.meta <- function(mua, mub, sa, sb, rho, n = NULL, alpha=0.05){
  pi <- 3.14159265359
  tt<-seq(0, 2*pi, length=1000)
  #CC<-sqrt(qchisq(1-alpha, 2))
  #CC<-sqrt(2*qf(1-alpha, 2, n-2))
  CC <- qnorm(1-alpha/2)
  ua <<- mua + sa*CC*cos(tt)
  ub <<- mub + sb*CC*cos(tt+acos(rho))
  return(data.frame(x = ua,y = ub))
}

# ellipse.plot.pade0: old code to draw ellipse
# source: https://stackoverflow.com/questions/41820683/how-to-plot-ellipse-given-a-general-equation-in-r
ellipse.plot.pade0 <- function(coef, ebar){
  
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


ellipse.plot.pade <- function(coef, ebar){
  
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