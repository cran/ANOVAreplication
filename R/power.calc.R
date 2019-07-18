#function to calculate power post-hoc with Ha: mu1 = mu.. = muJ
power.calc <- function(n.r,posterior,g.m,p.sd,
                       statistic,Amat=0L,exact=0L,difmin=0L,effectsize=FALSE,
                       alpha=.05){
  p <- ncol(posterior)-1    #number of groups
  lFps <- dim(posterior)[1] #length future Fps

  #null distribution = F scores for original dataset
  power.H0 <- prior.predictive.check(n=n.r,posterior=posterior,obs=FALSE,
                                     statistic=statistic,Amat=Amat,exact=exact,difmin=difmin,effectsize=effectsize)
  #rejection value is max 5% of H0
  Fps.power.H0 <- power.H0$F_sim
  rej.value <- quantile(Fps.power.H0, 1-alpha)

  #alternative distr. = F scores when all group means are equal (g.m = general mean original)
  #the SD of the y-data is that of the original dataset (p.sd)
  power.H1 <- prior.predictive.check(n=n.r,posterior=cbind(matrix(g.m,nrow=lFps,ncol=p,byrow=TRUE),rep(p.sd,lFps)),
                                     obs=FALSE,statistic=statistic,Amat=Amat,exact=exact,difmin=difmin,effectsize=effectsize)
  Fps.power.H1 <- power.H1$F_sim

  #power, proportion F's more extreme than rej value based on null distribution (F's original)
  power.out <- sum(Fps.power.H1>rej.value)/length(Fps.power.H1)

  hist(Fps.power.H0,freq=FALSE,col=rgb(1,0,0,1/4),border=rgb(1,0,0,1/2),
       breaks=c(seq(0,max(Fps.power.H0)+10)),main="",xlab=expression(bar(F))) #null with true effect
  hist(Fps.power.H1,freq=FALSE,col=rgb(0,0,1,1/4),border=rgb(0,0,1,1/2),
       breaks=c(seq(0,max(Fps.power.H1)+10)),add=TRUE) #H1, means equal
  abline(v=rej.value,col=rgb(1,0,1,1/2),lwd=2) #power is blue at the right side of this line

  return.info <- list(power=power.out,rejection.value=rej.value,F.H0=Fps.power.H0,F.H1=Fps.power.H1)
  return(return.info)
}
