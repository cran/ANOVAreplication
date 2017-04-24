#prior predictive check ####
#function with input N: replication study sample size, n: subgroup sample sizes,
#posterior: samples from posterior original data,
#statistic: Fbar inequality (ineq), Fbar exact valuees (exact), specific differences (dif),
#obs: observed replication statistic provided TRUE for replication but FALSE for power check,
#F_n: actual observed F observed in replication.
#repdata: whether datasets are already generated in a previous run of the check, default = FALSE
#Amat: matrix with inequality constraints for statistic = "ineq" and "dif".
#exact: vector with exact mean values for statistic = "exact"
#difmin: vector with hypothesized differences for statistic = "dif"
prior.predictive.check <- function(n,posterior,statistic,obs=TRUE,F_n,
                                   repdata=FALSE,Amat=0L,exact=0L,difmin=0L,effectsize=FALSE,
                                   seed=0){
  if(seed!=0){
    set.seed(seed)}
  N=sum(n)
  p <- length(n)         #number of groups
  it <- dim(posterior)[1]
  if(repdata==FALSE){ #if p(y) is not already generated
    y <- array(NA, dim=c(N,1,it))
    x <- array(NA, dim=c(N,1,it))
    data.rep <- list()
    for (l in 1:it){
      #create data from normal prior based on sampled values
      for (i in 1:p){
        x[(sum(n[1:i])+1-n[i]):sum(n[1:i]),1,l] <- c(rep(i,n[i]))
        y[(sum(n[1:i])+1-n[i]):sum(n[1:i]),1,l] <- c(rnorm(n=n[i],mean=posterior[l,i],sd=posterior[l,p+1]))}
      data.rep[[l]] <- as.data.frame(cbind(y[,1,l],x[,1,l]))}
  }else{data.rep=data.rep} #else, use existing p(y)
  #on each dataset, apply the function of interest and store results (Fbar & df-error) in Fps
  if(statistic=="ineq"){
    Fps <- unlist(lapply(X=data.rep, FUN = function(x){Fbar.ineq(x,Amat=Amat)}))}
  if(statistic=="exact"){
    Fps <- unlist(lapply(X=data.rep, FUN = function(x){Fbar.exact(x,exact=exact)}))}
  if(statistic=="dif"){
    Fps <- unlist(lapply(X=data.rep, FUN = function(x){Fbar.dif(x,Amat=Amat,difmin=difmin)}))}
  if(statistic=="dif"&effectsize==TRUE){
    Fps <- unlist(lapply(X=data.rep, FUN = function(x){Fbar.dif(x,Amat=Amat,difmin=difmin, effectsize=TRUE)}))}
  if(obs==TRUE){
    #calculate proportion of replicated F more extreme than observed = p-value, one-sided
    ppp <- sum(Fps>=F_n)/it

    #make plot
    par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
    hist(Fps,main="",xlab="",freq=TRUE,las=1)
    abline(v=F_n,col="red")

    result <- list("sumFdist"=summary(Fps),"ppp"=round(ppp,4),"F_sim"=Fps)
    return(result)}

  if(obs==FALSE){
    result <- list("F_sim"=Fps)
    return(result)
  }
}
