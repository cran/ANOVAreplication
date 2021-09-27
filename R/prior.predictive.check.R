#prior predictive check ####
#function with input N: replication study sample size, n: subgroup sample sizes,
#posterior: samples from posterior original data,
#statistic: Fbar inequality (ineq), Fbar exact valuees (exact), specific differences (dif),
#obs: observed replication statistic provided TRUE for replication but FALSE for power check,
#F_n: actual observed F observed in replication.
#Amat: matrix with inequality constraints for statistic = "ineq" and "dif".
#exact: vector with exact mean values for statistic = "exact"
#difmin: vector with hypothesized differences for statistic = "dif"
prior.predictive.check <- function(n,posterior,statistic,obs=TRUE,F_n,
                                   Amat=0L,exact=0L,difmin=0L,effectsize=FALSE,
                                   seed=0){
  warning("The constraints should be in line with the original data. Test this by computing F-bar for the original data. The resulting F-bar should equal 0. ", call. = FALSE)

  if(seed!=0){
    set.seed(seed)}

  lapply_pb <- function(X, FUN, ...)
  {
    env <- environment()
    pb_Total <- length(X)
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

    # wrapper around FUN
    wrapper <- function(...){
      curVal <- get("counter", envir = env)
      assign("counter", curVal +1 ,envir=env)
      setTxtProgressBar(get("pb", envir=env), curVal +1)
      FUN(...)
    }
    res <- lapply(X, wrapper, ...)
    close(pb)
    res
  }

  N=sum(n)
  p <- length(n)         #number of groups
  it <- dim(posterior)[1]

    y <- array(NA, dim=c(N,1,it))
    x <- array(NA, dim=c(N,1,it))
    data.rep <- list()
    for (l in 1:it){
      #create data from normal prior based on sampled values
      for (i in 1:p){
        x[(sum(n[1:i])+1-n[i]):sum(n[1:i]),1,l] <- c(rep(i,n[i]))
        y[(sum(n[1:i])+1-n[i]):sum(n[1:i]),1,l] <- c(rnorm(n=n[i],mean=posterior[l,i],sd=posterior[l,p+1]))}
      data.rep[[l]] <- as.data.frame(cbind(y[,1,l],x[,1,l]))}

  #on each dataset, apply the function of interest and store results (Fbar & df-error) in Fps
  if(statistic=="ineq"){
    Fps <- unlist(lapply_pb(X=data.rep, FUN = function(x){Fbar.ineq(x,Amat=Amat)}))}
  if(statistic=="exact"){
    Fps <- unlist(lapply_pb(X=data.rep, FUN = function(x){Fbar.exact(x,exact=exact)}))}
  if(statistic=="dif"&effectsize==FALSE){
    Fps <- unlist(lapply_pb(X=data.rep, FUN = function(x){Fbar.dif(x,Amat=Amat,difmin=difmin)}))}
  if(statistic=="dif"&effectsize==TRUE){
    Fps <- unlist(lapply_pb(X=data.rep, FUN = function(x){Fbar.dif(x,Amat=Amat,difmin=difmin, effectsize=TRUE)}))}
  if(obs==TRUE){
    #calculate proportion of replicated F more extreme than observed = p-value, one-sided
    ppp <- sum(Fps>=F_n)/it

    #make plot
    #par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
    #hist(Fps,main="",xlab="",freq=TRUE,las=1)
    #abline(v=F_n,col="red")

    #make plot 2
    #save and restore user par
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    #set par for output plots
    par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
    Mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]}

    if(sum(Fps==Mode(Fps))/length(Fps)>=.10){
      h <- hist(Fps[-which(Fps==Mode(Fps))],ylim=c(0,sum(Fps==Mode(Fps))),xlab=expression(italic(bar(F)[bold(y)][r])),ylab="Frequency",main="",breaks=20)
      hist(Fps[-which(Fps==Mode(Fps))],ylim=c(0,max(sum(Fps==Mode(Fps)),h$counts[1])),
           xlim=c(0,max(max(Fps),F_n)),
           xlab=expression(italic(bar(F)[bold(y)])),ylab="Frequency",main="",breaks=20)
      segments(x0=Mode(Fps),y0=0,x1=Mode(Fps),y1=sum(Fps==Mode(Fps)),col="black",lwd=5)
      abline(v=F_n,col="red")
    }else{
      hist(Fps,freq=TRUE,breaks=seq(0,max(Fps),length.out=40),xlim=c(0,max(max(Fps),F_n)),
           xlab=expression(italic(bar(F)[bold(y)])),ylab="Frequency",main="")
      abline(v=F_n,col="red")}

    result <- list("sumFdist"=summary(Fps),"ppp"=round(ppp,4),"F_sim"=Fps)
    return(result)}

  if(obs==FALSE){
    result <- list("F_sim"=Fps)
    return(result)
  }
}
