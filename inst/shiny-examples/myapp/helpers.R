#https://osf.io/fgjvw/ = all replications
#https://osf.io/pz0my/ = specific project Replication of Monin, Sawyer, & Marquez (2008, JPSP 95(1), Exp. 4)
#analysis: https://darrenjw.wordpress.com/2014/12/22/one-way-anova-with-fixed-and-random-effects-from-a-bayesian-perspective/
#repproject <- read.csv("rpp_data.csv")
#repMonin <- repproject[repproject$Study.Num==43,]

#function to generate exact data ####
rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }

#F-bar functions ####
library(quadprog)

#Fbar function, input = data.frame, Amat is a matrix with inequality constraints,
#each inequality is summarized within a line that has k indices. For example u1 < u2, u3 = rbind(c(-1,1,0),c(-1,0,1)),
#Fresult contains the resulting Fbar, and the RSS-df
Fbar.ineq <- function(data,Amat){
  names(data) <- c("V1","V2")
  fit.lm <- lm(V1~as.factor(V2)-1,data)                   #lm, no intercept (!)
  mfit <- fit.lm$model                                    
  Y <- model.response(mfit)                               #the data stored in Y
  X <- model.matrix(fit.lm)[,,drop = FALSE]               #dummies
  s2 <- summary(fit.lm)$sigma^2                           #Residual Standard Error squared
  df.error <- summary(fit.lm)$fstatistic[[3]]             #error df
  XX <- crossprod(X); Xy <- t(X) %*% Y
  out.h0 <- solve.QP(Dmat = XX, dvec = Xy, Amat = t(Amat))
  RSS.h0 <- sum((Y - (X %*% out.h0$solution))^2)
  RSS.ha <- sum((Y - (X %*% out.h0$unconstrained.solution))^2)
  #hypothesis test Type B 
  Fresult <<- (RSS.h0 - RSS.ha)/s2
  return(Fresult)
}

#Fbar for exact values. For conceptual replication, use standardized dependent, and standardized exact values
#data = data.frame, exact=the exact values for each of the parameters 
#possible extension: exact value for some parameters
Fbar.exact <- function(data,exact){
  names(data) <- c("V1","V2")
  fit.lm <- lm(V1~as.factor(V2)-1,data)                   #lm, no intercept (!)
  mfit <- fit.lm$model                                    #standard linear model
  Y <- model.response(mfit)                               #the data stored in Y
  X <- model.matrix(fit.lm)[,,drop = FALSE]               #dummies, intercept + K-1 groups
  s2 <- summary(fit.lm)$sigma^2                           #Residual Standard Error squared
  df.error <- summary(fit.lm)$fstatistic[[3]]             #error df
  RSS.h0 <- sum((Y - (X %*% exact))^2)
  RSS.ha <- sum((Y - (X %*% fit.lm$coefficients))^2)
  #hypothesis test Type B 
  Fresult <<- (RSS.h0 - RSS.ha)/s2
  return(Fresult)
}

#Fbar for hypothesis with minimal (directional) difference between means. data = data.frame
#Amat is the matrix with inequality constraints, difmin is the minimal difference per Amat row in vector form
Fbar.dif <- function(data,Amat, difmin, effectsize=FALSE){
  names(data) <- c("V1","V2")
  fit.lm <- lm(V1~as.factor(V2)-1,data)                   #lm, no intercept (!)
  mfit <- fit.lm$model                                    #standard linear model
  Y <- model.response(mfit)                               #the data stored in Y
  X <- model.matrix(fit.lm)[,,drop = FALSE]               #dummies, intercept + K-1 groups
  s2 <- summary(fit.lm)$sigma^2                           #Residual Standard Error squared
  df.error <- summary(fit.lm)$fstatistic[[3]]             #error df
  XX <- crossprod(X); Xy <- t(X) %*% Y
  if(effectsize==TRUE){
    a <- apply(Amat, MARGIN=1, FUN= function(x){which(x==-1)})
    b <- apply(Amat, MARGIN=1, FUN= function(x){which(x==1)})
    n.r <- as.numeric(table(data[,2]))
    s <- unlist(lapply(1:dim(Amat)[1], FUN=function(x,i){
      sqrt(((n.r[a[i]]-1)*var(x$V1[which(x$V2==a[i])])+(n.r[b[i]]-1)*var(x$V1[which(x$V2==b[i])]))/
             (n.r[a[i]]+n.r[b[i]]-2))},x=as.list(data)))    
    out.h0 <- solve.QP(Dmat = XX, dvec = Xy, Amat = t(Amat),bvec=difmin*s)
  }else{
    out.h0 <- solve.QP(Dmat = XX, dvec = Xy, Amat = t(Amat),bvec=difmin)
  }  
  #print(out.h0$solution) #print(out.h0$unconstrained.solution)
  RSS.h0 <- sum((Y - (X %*% out.h0$solution))^2)
  RSS.ha <- sum((Y - (X %*% out.h0$unconstrained.solution))^2)
  #hypothesis test Type B 
  Fresult <<- (RSS.h0 - RSS.ha)/s2
  return(Fresult)
}


#prior predicitve check ####
#function with input N: replication study sample size, n: subgroup sample sizes,
#posterior: samples from posterior original data,
#statistic: F classic (null), Fbar inequality (ineq), Fbar exact valuees (exact), specific differences (dif),
#obs: observed replication statistic provided TRUE for replication but FALSE for power check,
#F_obs: actual observed F observed in replication.
#repdata: whether datasets are already generated in a previous run of the check, default = FALSE
#Amat: matrix with inequality constraints for statistic = "ineq" and "dif".
#exact: vector with exact mean values for statistic = "exact"
#difmin: vector with hypothesized differences for statistic = "dif"
prior.predictive.check <- function(n,posterior,statistic,obs=TRUE,F_obs,
                                   repdata=FALSE,Amat=0L,exact=0L,difmin=0L,effectsize=FALSE,
                                   seed=0){
  if(seed!=0){
    set.seed(seed=seed)}
  N=sum(n)
  p <- length(n)         #number of groups
  it <- dim(posterior)[1]
  if(repdata==FALSE){ #if p(y) is not already generated
    y <- array(NA, dim=c(N,1,it))
    x <- array(NA, dim=c(N,1,it))
    data.rep <- list()
    # Create a Progress object
    progress <- shiny::Progress$new()
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())

    progress$set(message = "Simulating data", value = 0)
    for (l in 1:it){
      #create data from normal prior based on sampled values
      if(l %% 1000==0) {progress$inc(1/it, detail = paste("Sample", l))}
      for (i in 1:p){
        x[(sum(n[1:i])+1-n[i]):sum(n[1:i]),1,l] <- c(rep(i,n[i]))
        y[(sum(n[1:i])+1-n[i]):sum(n[1:i]),1,l] <- c(rnorm(n=n[i],mean=posterior[l,i],sd=posterior[l,p+1]))}
      data.rep[[l]] <- as.data.frame(cbind(y[,1,l],x[,1,l]))}
  }else{data.rep=data.rep} #else, use existing p(y)
  #on each dataset, apply the function of interest and store results (Fbar & df-error) in Fps
  if(statistic=="ineq"){
    Fps <<- matrix(unlist(
      lapply(X=data.rep, FUN = function(x){Fbar.ineq(x,Amat=Amat)})),ncol=1,byrow=TRUE)}
  if(statistic=="exact"){
    Fps <<- matrix(unlist(
      lapply(X=data.rep, FUN = function(x){Fbar.exact(x,exact=exact)})),ncol=1,byrow=TRUE)}
  if(statistic=="dif"&effectsize==FALSE){
    Fps <<- matrix(unlist(
      lapply(X=data.rep, FUN = function(x){Fbar.dif(x,Amat=Amat,difmin=difmin)})),ncol=1,byrow=TRUE)}
  if(statistic=="dif"&effectsize==TRUE){
    Fps <<- matrix(unlist(
      lapply(X=data.rep, FUN = function(x){Fbar.dif(x,Amat=Amat,difmin=difmin,effectsize=TRUE)})),ncol=1,byrow=TRUE)}
  if(obs==TRUE){
    #calculate proportion of replicated F more extreme than observed = p-value, one-sided
    ppp <<- sum(Fps>=F_obs)/it
    Fps.summary <- as.matrix(summary(Fps),ncol=1)

    result <- list("distribution F-bar given original data"=Fps.summary,"F-bar replication data"=F_obs,"prior predictive p-value"=round(ppp,4))
    return(result)}
}

#power ####
power.ppp <- function(start_n,powtarget=.825,powmargin=.025,posterior,g.m,statistic,
                      Amat=0L,exact=0L,difmin=0L,effectsize=FALSE,
                      nmax=600,alpha=.05,itmax=10){
  it=0; power.out=0
  Npower.l <<- matrix(NA,ncol=2,nrow=itmax)
  p <- ncol(posterior)-1    #number of groups
  lFps <- dim(posterior)[1] #length future Fps
  Npower = p*start_n

  exit=FALSE
  while(exit==FALSE&&it<itmax&&
        (power.out<(powtarget-powmargin)|power.out>(powtarget+powmargin))){
    it=it+1

    if(it>itmax){ #when max iteration, stop before next loop (exit = TRUE)
      stop <- "The maximum number of iterations has been reached"
      exit=TRUE
      }
    if(Npower>=nmax){ #when max sample size is reached, stop before next loop (exit = TRUE)
      Npower = nmax
      stop <- "The total sample size was reset to its maximum (default = 600) to limit computational time"
      exit=TRUE
      }

    if(Npower<20){    #min sample size, avoiding negative sample sizes
      print("The total sample size was reset to its minimum of 20")
      Npower = 20}

    nF <- ceiling(c(rep(Npower/p,p)))   #subsample sizes, equal over groups
    Npower <- sum(nF)                   #actual sample size after equal groups

    #null distribution = F scores for original dataset
    prior.predictive.check(n=nF,posterior=posterior,obs=FALSE,
                           statistic=statistic,Amat=Amat,exact=exact,difmin=difmin,effectsize=effectsize)
    #rejection value is max 5% of H0
    Fps.power.H0 <- Fps
    rej.value <<- quantile(Fps.power.H0, 1-alpha)

    #alternative distr. = F scores when all group means are equal (value = general mean original)
    #the SE of posterior/prior means is derived from the original posterior
    #the SD of the y-data is that of the original posterior
    posteriormeans.A <- matrix(NA,nrow=lFps,ncol=p)
    for (i in 1:p){posteriormeans.A[,i] <- rnorm(lFps,mean=g.m,sd=apply(posterior[,1:p],2,sd))}
    prior.predictive.check(n=nF,posterior=cbind(posteriormeans.A,posterior[,p+1]),
                           obs=FALSE,statistic=statistic,Amat=Amat,exact=exact,difmin=difmin,effectsize=effectsize)

    #power, proportion F's more extreme than rej value based on null distribution (F's original)
    Fps.power.H1 <- Fps
    power.out <- sum(Fps.power.H1>rej.value)/length(Fps.power.H1)
    Npower.l[it,] <<- c(nF[1],round(power.out,2)) #store sample size per group and power

    #calculations sample size for next loop
    #a1) after two iterations, calculate necessary N with linear regression
    #a2) if slope is negative, use method b
    #b) multiply sample size with x, where x is the ratio of current power and target power
    if(it>1){
      #exit conditions that can be observed after second iteration
      if(round(Npower.l[it,2],2)==round(Npower.l[it-1,2],2)){
        stop <- "No change in power. Potentially a maximum power level is reached. Limitations to power for simple order restrictions can be further explored with the R-package complexity"
        exit=TRUE}
      if(Npower.l[it,1]==Npower.l[it-1,1]){
        stop <- "No sample size variation"
        exit=TRUE }}

      if(it>1&&power.out<(powtarget-powmargin)){

      if(it==3){          #quadratic regression
        pow.coef <- lm(Npower.l[,2]~Npower.l[,1]+I(Npower.l[,1]^2))$coefficients #power coefficients
        Npower <- Re(suppressWarnings(polyroot(c(pow.coef[1]-powtarget,pow.coef[2:3])))[1])*p}else{    #linear regression
            ipu <- lm(Npower.l[,2]~Npower.l[,1])$coefficients[2] #power increase per unit
            if(ipu>0){                                           #power should only increase with >n
              itarget <- powtarget-power.out
              Npower = (nF[1]+as.numeric(itarget/ipu))*p}else{
                x <- powtarget/power.out                   #ratio current power, target power
                Npower <- Npower*x}
          }
    }else{
      x <- powtarget/power.out                   #ratio current power, target power
      Npower <- Npower*x}
  }
  colnames(Npower.l) <- c("n per group","Power")
  Fps.power.H0 <<- Fps.power.H0
  Fps.power.H1 <<- Fps.power.H1

  if (it == itmax){stop <- "The maximum number of iterations has been reached."}
  if (power.out>(powtarget-powmargin)&power.out<(powtarget+powmargin)){
    stop <- "The target power level has been reached."}

  return.info <- list(stop,Npower.l[1:it,])
  return(return.info)}
