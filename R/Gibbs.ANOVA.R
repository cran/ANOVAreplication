
Gibbs.ANOVA <- function(data,it=5000,burnin=500){
  #R program for Gibbs sampling from conditionals

  I=it+burnin                                      #number of iterations
  #read only observations with complete information
  x=data$g
  y=data$y
  N=length(y)
  G=length(unique(data$g))

  #other file input
  fit.lm <- lm(y~as.factor(x)-1)                   #lm, no intercept (!)
  x <- model.matrix(fit.lm)[,,drop = FALSE]

  #establish parameter vectors and constant quantities
  s1=matrix(1,I); b1=matrix(0,I,G)
  s2=matrix(1,I); b2=matrix(0,I,G)
  xtxi=solve(t(x)%*%x)
  pars=coefficients(lm(y~x-1))
  #Gibbs sampling begins
  for(i in 2:I){
    #simulate beta from its multivariate normal conditional
    b1[i,]=pars+t(rnorm(G,mean=0,sd=1))%*%chol((s1[i-1]^2)*xtxi) #choleski decomposition
    #simulate sigma from its inverse gamma distribution
    s1[i]=sqrt(1/rgamma(1,N/2,.5*t(y-x%*%(b1[i,]))%*%(y-x%*%(b1[i,]))))
  }
  for(i in 2:I){
    #simulate beta from its multivariate normal conditional
    b2[i,]=pars+t(rnorm(G,mean=0,sd=1))%*%chol((s2[i-1]^2)*xtxi) #choleski decomposition
    #simulate sigma from its inverse gamma distribution
    s2[i]=sqrt(1/rgamma(1,N/2,.5*t(y-x%*%(b2[i,]))%*%(y-x%*%(b2[i,]))))
  }
  ## plot posterior density for parameters of interest with indications for credible interval
  par1=cbind(b1,s1)[-c(1:burnin),]
  par2=cbind(b2,s2)[-c(1:burnin),]
  par=list(par1,par2)
  output_m <<- rbind(par1,par2)
  colnames(output_m) <- paste("Mean",1:(G+1))
  colnames(output_m)[G+1] <- "SD"

  par(mfrow = c(G+1, 2), mar = rep(2, 4))

  for (j in 1:(G+1)) {

    hist(output_m[,j], breaks="Scott", main=colnames(par)[j], xlab="Value")
    abline(v=quantile(output_m[,j], probs=c(0.025, 0.975)), col="blue")

    plot.ts(cbind(par1[,j],par2[,j]),col=c("blue","green"), plot.type="single",ylab="")

  }

results <- cbind(round(apply(output_m,2,mean),4),round(apply(output_m,2,median),4),
                 round(apply(output_m,2,sd),4))
colnames(results) <- c("mean","median","sd")
  return(results)
}
