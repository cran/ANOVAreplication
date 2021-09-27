sample.size.calc <- function(start_n=20,itmax=10,nmax=600,powtarget=.825,powmargin=.025,
                             posterior,g.m,p.sd,
                             statistic,Amat=0L,exact=0L,difmin=0L,effectsize=FALSE,
                             alpha=.05,printit=TRUE){
  it=0; power.out=0
  Npower.l <- matrix(NA,ncol=2,nrow=itmax)
  colnames(Npower.l) <- c("n per group","Power")
  p <- ncol(posterior)-1    #number of groups
  lFps <- dim(posterior)[1] #length future Fps
  Npower = p*start_n

  exit=FALSE
  while(exit==FALSE&&it<itmax&&
        (power.out<(powtarget-powmargin)|power.out>(powtarget+powmargin))){
    it=it+1
    if(printit==TRUE){print(it)}
    if(it>itmax){ #when max iteration, stop before next loop (exit = TRUE)
      warning("The maximum number of iterations has been reached")
      exit=TRUE}
    if(Npower>=nmax){ #when max sample size is reached, stop before next loop (exit = TRUE)
      Npower = nmax
      warning("The total sample size was reset to its maximum (default = 600) to limit computational time")
      exit=TRUE}
    if(Npower<20){    #min sample size, avoiding negative sample sizes
      Npower = 20
      print(Npower.l)
      warning("The total sample size was reset to its minimum of 20")}
    nF <- ceiling(c(rep(Npower/p,p)))   #subsample sizes, equal over groups
    Npower <- sum(nF)                   #actual sample size after equal groups

    #null distribution = F scores for original dataset

    H0 <- suppressWarnings(prior.predictive.check(n=nF,posterior=posterior,obs=FALSE,
                                 statistic=statistic,Amat=Amat,exact=exact,difmin=difmin,effectsize=effectsize))
    #rejection value is max 5% of H0
    F_sim.H0 <- H0$F_sim
    rej.value <- quantile(F_sim.H0, 1-alpha)

    #alternative distr. = F scores when all group means are equal (value = general mean original)
    #the pooled SD of the y-data is that of the original posterior
    Ha <- suppressWarnings(prior.predictive.check(n=nF,posterior=cbind(matrix(g.m,nrow=lFps,ncol=p,byrow=TRUE),rep(p.sd,lFps)),
                                 obs=FALSE,statistic=statistic,Amat=Amat,exact=exact,difmin=difmin,effectsize=effectsize))

    #power, proportion F's more extreme than rej value based on null distribution (F's original)
    F_sim.Ha <- Ha$F_sim
    power.out <- sum(F_sim.Ha>rej.value)/length(F_sim.Ha)
    Npower.l[it,] <- c(nF[1],round(power.out,2)) #store sample size per group and power

    #calculations sample size for next loop
    #a1) after two iterations, calculate necessary N with linear regression
    #a2) if slope is negative, use method b
    #b) multiply sample size with x, where x is the ratio of current power and target power
    if(it>1&&power.out<(powtarget-powmargin)){
      if(round(Npower.l[it,2],2)==round(Npower.l[it-1,2],2)){exit=TRUE;
      warning("No change in power. Potentially a maximum power level is reached. Limitations to power can be further explored with the complexity function")}
      if(Npower.l[it,1]==Npower.l[it-1,1]){exit=TRUE; warning("No sample size variation")}

      if(it==3){
        pow.coef <- lm(Npower.l[,2]~Npower.l[,1]+I(Npower.l[,1]^2))$coefficients #power coefficients
        Npower <- Re(polyroot(c(pow.coef[1]-powtarget,pow.coef[2:3]))[1])*p}else{

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
  hist(F_sim.H0,freq=FALSE,col=rgb(1,0,0,1/4),border=rgb(1,0,0,1/2),main="",xlab=expression(bar(F))) #null with true effect
  hist(F_sim.Ha,freq=FALSE,col=rgb(0,0,1,1/4),border=rgb(0,0,1,1/2),add=TRUE) #H1, means equal
  abline(v=rej.value,col=rgb(1,0,1,1/2),lwd=2)
  return(list(sampcalc=Npower.l[1:it,]))
  }
