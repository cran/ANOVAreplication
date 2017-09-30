#F-bar functions ####

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
  XX <- crossprod(X); Xy <- t(X) %*% Y
  out.h0 <- solve.QP(Dmat = XX, dvec = Xy, Amat = t(Amat))
  RSS.h0 <- sum((Y - (X %*% out.h0$solution))^2)
  RSS.ha <- sum((Y - (X %*% out.h0$unconstrained.solution))^2)
  #hypothesis test Type B
  Fresult <- c((RSS.h0 - RSS.ha)/s2)
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
  RSS.ha <- sum(resid(fit.lm)^2)
  #hypothesis test Type B
  Fresult <- (RSS.h0 - RSS.ha)/s2
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
  Fresult <- (RSS.h0 - RSS.ha)/s2
  return(Fresult)
}

