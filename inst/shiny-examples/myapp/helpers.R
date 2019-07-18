#https://osf.io/fgjvw/ = all replications
#https://osf.io/pz0my/ = specific project Replication of Monin, Sawyer, & Marquez (2008, JPSP 95(1), Exp. 4)
#analysis: https://darrenjw.wordpress.com/2014/12/22/one-way-anova-with-fixed-and-random-effects-from-a-bayesian-perspective/
#repproject <- read.csv("rpp_data.csv")
#repMonin <- repproject[repproject$Study.Num==43,]

#function to generate exact data ####
generate.data <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }

pooled.sd <- function(data){
  p <- length(table(data$g))
  n.g <- table(data$g)
  sd.g <- aggregate(data$y,by=list(data$g),sd)[,2]
  sqrt(sum((n.g - 1) * sd.g^2)/(sum(n.g) - p))
}

#F-bar functions ####
library(quadprog)

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
  Fresult <- (RSS.h0 - RSS.ha)/s2
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


#prior predicitve check ####
#function with input N: replication study sample size, n: subgroup sample sizes,
#posterior: samples from posterior original data,
#statistic: F classic (null), Fbar inequality (ineq), Fbar exact valuees (exact), specific differences (dif),
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
    set.seed(seed=seed)}
  N=sum(n)
  p <- length(n)         #number of groups
  it <- dim(posterior)[1]
  if(repdata==FALSE){ #if p(y) is not already generated
    y <- array(NA, dim=c(N,1,it))
    x <- array(NA, dim=c(N,1,it))
    data.rep <- list()
    # Create a Progress object
    progress <- shiny::Progress$new(min=0,max=it)
    # Make sure it closes when we exit this reactive, even if there's an error
    on.exit(progress$close())

    progress$set(message = "Simulating data", value = 0)
    for (l in 1:it){
      #create data from normal prior based on sampled values
      if(l %% 1000==0) {progress$inc(it, detail = paste("Sample", l))}
      for (i in 1:p){
        x[(sum(n[1:i])+1-n[i]):sum(n[1:i]),1,l] <- c(rep(i,n[i]))
        y[(sum(n[1:i])+1-n[i]):sum(n[1:i]),1,l] <- c(rnorm(n=n[i],mean=posterior[l,i],sd=posterior[l,p+1]))}
      data.rep[[l]] <- as.data.frame(cbind(y[,1,l],x[,1,l]))}
  }else{data.rep=data.rep} #else, use existing p(y)
  #on each dataset, apply the function of interest and store results (Fbar & df-error) in Fps

  lapply_pb <- function(X, FUN, ...)
  {
    env <- environment()
    pb_Total <- length(X)
    counter <- 0
    pb <- Progress$new(min=0, max = pb_Total)

    # wrapper around FUN
    wrapper <- function(...){
      curVal <- get("counter", envir = env)
      pb$set(value = curVal +1)
      if(curVal %% 100==0) {
        pb$inc(curVal, message = "Calculating", detail = paste("Fbar", curVal))
      }
      assign("counter", curVal +1 ,envir=env)
      #setTxtProgressBar(get("pb", envir=env), curVal +1)
      FUN(...)
    }
    res <- lapply(X, wrapper, ...)
    on.exit(pb$close())
    res
  }

  if(statistic=="ineq"){
    Fps <<- matrix(unlist(
      lapply_pb(X=data.rep, FUN = function(x){Fbar.ineq(x,Amat=Amat)})),ncol=1,byrow=TRUE)}
  if(statistic=="exact"){
    Fps <<- matrix(unlist(
      lapply_pb(X=data.rep, FUN = function(x){Fbar.exact(x,exact=exact)})),ncol=1,byrow=TRUE)}
  if(statistic=="dif"&effectsize==FALSE){
    Fps <<- matrix(unlist(
      lapply_pb(X=data.rep, FUN = function(x){Fbar.dif(x,Amat=Amat,difmin=difmin)})),ncol=1,byrow=TRUE)}
  if(statistic=="dif"&effectsize==TRUE){
    Fps <<- matrix(unlist(
      lapply_pb(X=data.rep, FUN = function(x){Fbar.dif(x,Amat=Amat,difmin=difmin,effectsize=TRUE)})),ncol=1,byrow=TRUE)}
  if(obs==TRUE){
    #calculate proportion of replicated F more extreme than observed = p-value, one-sided
    ppp <<- sum(Fps>=F_n)/it
    Fps.summary <- as.matrix(summary(Fps),ncol=1)

    result <- list("distribution F-bar given original data"=Fps.summary,"F-bar replication data"=F_n,"prior predictive p-value"=round(ppp,4),
                   "F_sim"=Fps)
    return(result)}
}

#power ####
#function to calculate power post-hoc with Ha: mu1 = mu.. = muJ
power.calc <- function(n.r,posterior,g.m,p.sd,
                       statistic,Amat=0L,exact=0L,difmin=0L,effectsize=FALSE,
                       alpha=.05){
  p <- ncol(posterior)-1    #number of groups
  lFps <- dim(posterior)[1] #length future Fps

  #null distribution = F scores for original dataset
  power.H0 <<- prior.predictive.check(n=n.r,posterior=posterior,obs=FALSE,
                                      statistic=statistic,Amat=Amat,exact=exact,difmin=difmin,effectsize=effectsize)
  #rejection value is max 5% of H0
  Fps.power.H0 <- Fps
  rej.value <- quantile(Fps.power.H0, 1-alpha)

  #alternative distr. = F scores when all group means are equal (value = general mean original)
  #the SD of the y-data is that of the original posterior
  power.H1 <<- prior.predictive.check(n=n.r,posterior=cbind(matrix(g.m,nrow=lFps,ncol=p,byrow = TRUE),rep(p.sd,lFps)),
                                      obs=FALSE,statistic=statistic,Amat=Amat,exact=exact,difmin=difmin,effectsize=effectsize)

  #power, proportion F's more extreme than rej value based on null distribution (F's original)
  Fps.power.H1 <- Fps

  #power, proportion F's more extreme than rej value based on null distribution (F's original)
  power.out <- sum(Fps.power.H1>rej.value)/length(Fps.power.H1)

  Fps.power.H0 <<- Fps.power.H0
  Fps.power.H1 <<- Fps.power.H1
  rej.value <<- rej.value

  return.info <- list(power=power.out,rejection.value=rej.value)
  return(return.info)
}

#sample size calculator ####
sample.size.calc <- function(start_n,itmax=10,nmax=600,powtarget=.825,powmargin=.025,
                             posterior,g.m,p.sd,
                             statistic,Amat=0L,exact=0L,difmin=0L,effectsize=FALSE,
                             alpha=.05){
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

    n.r <- ceiling(c(rep(Npower/p,p)))   #subsample sizes, equal over groups
    Npower <- sum(n.r)                   #actual sample size after equal groups

    #null distribution = F scores for original dataset
    power.H0 <<- prior.predictive.check(n=n.r,posterior=posterior,obs=FALSE,
                                 statistic=statistic,Amat=Amat,exact=exact,difmin=difmin,effectsize=effectsize)
    #rejection value is max 5% of H0
    Fps.power.H0 <<- Fps
    rej.value <<- quantile(Fps.power.H0, 1-alpha)

    #alternative distr. = F scores when all group means are equal (value = general mean original)
    #the SD of the y-data is that of the original posterior
    power.Ha <<- prior.predictive.check(n=n.r,posterior=cbind(matrix(g.m,nrow=lFps,ncol=p,byrow = TRUE),rep(p.sd,lFps)),
                                 obs=FALSE,statistic=statistic,Amat=Amat,exact=exact,difmin=difmin,effectsize=effectsize)

    #power, proportion F's more extreme than rej value based on null distribution (F's original)
    Fps.power.H1 <<- Fps
    power.out <- sum(Fps.power.H1>rej.value)/length(Fps.power.H1)
    Npower.l[it,] <<- c(n.r[1],round(power.out,2)) #store sample size per group and power

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
            Npower = (n.r[1]+as.numeric(itarget/ipu))*p}else{
              x <- powtarget/power.out                   #ratio current power, target power
              Npower <- Npower*x}
        }
    }else{
      x <- powtarget/power.out                   #ratio current power, target power
      Npower <- Npower*x}
  }
  Fps.power.H0 <<- Fps.power.H0
  Fps.power.H1 <<- Fps.power.H1

  if (it == itmax){stop <- "The maximum number of iterations has been reached."}
  if (power.out>(powtarget-powmargin)&power.out<(powtarget+powmargin)){
    stop <- "The target power level has been reached."}

  colnames(Npower.l) <- c("n per group","Power") #do not change location
  return.info <- list(stop,Npower.l[1:it,])
  return(return.info)}


#create matrices ####

#' Create (in)equality constraint matrices
#'
#' Parses a character string describing an informative hypothesis,
#' and returns (in)equality constraint matrices
#'
#' Informative hypotheses specified as a character string by "hyp" should
#' adhere to the following simple syntax: \itemize{
#' \item The hypothesis consists of a (series of) (in)equality
#' constraint(s). Every single (in)equality constraint is of the form "R1*mu1 +
#' R2*mu2+... = r", where capital Rs refer to numeric scaling constants, must
#' refer to the names of parameters in the model, and the lower case r refers
#' to a constant. Standard mathematical simplification rules apply; thus,
#' "R1*mu1 = R2*mu2" is equivalent to "R1*mu1 - R2*mu2 = 0".  \item Multiple
#' unrelated constraints within one hypothesis can be chained by "&". Thus,
#' "a=b&c=d" means that H1: a=b AND c=d.  \item Multiple related constraints
#' within one hypothesis can be chained by repeating the (in)equality operators
#' "=", "<", or ">". Thus, "a<b<c" means that H1: a < b AND b < c.  \item
#' Parameters can be grouped by placing them in a parenthesized, comma
#' separated list. Thus, "(a,b)>c" means that H1: a > c AND b > c.  Similarly,
#' "(a,b)>(c,d)" means that H1: a > c AND b > c AND b > c AND b > d.  }
#'
#' @aliases create_matrices create_matrices
#' @param varnames A character (vector of characters), containing names of
#' variables used in the hypotheses.  %Object of class \code{\link{lm}}, from
#' which the model parameters are extracted.
#' @param hyp A character string, containing a Bain hypothesis (see Details).
#' @return A pair of named matrices for every hypothesis specified in the
#' \code{hyp} argument; one matrix named ERr, specifying equality constraints,
#' and one matrix named IRr, specifying inequality constraints.
#' @author Caspar van Lissa
#' @keywords internal utilities
create_matrices <- function(varnames, hyp){
  if(is.null(varnames)) stop("Please input proper linear model object")
  hyp <- gsub("\\s", "", hyp)
  if(grepl("[><=]{2,}", hyp)) stop("Do not use combined comparison signs e.g., '>=' or '==', and use '&' to add a related constraint instead of ','")

  hyp_list <- strsplit(hyp, ";")[[1]]                         #mz
  hyp_list <- lapply(hyp_list, function(x){ strsplit(x, "&")[[1]]})
  hyp_list <- lapply(hyp_list, function(x){unlist(lapply(x, expand_compound_constraints))})
  hyp_list <- lapply(hyp_list, function(x){unlist(lapply(x, expand_parentheses))})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, flip_inequality)})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, constraint_to_equation)})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, order_terms)})

  hyp_list <- unlist(lapply(hyp_list, function(x){
    ERr <- x[grep("=", x)]
    IRr <- x[grep("[<>]", x)]
    if(length(ERr) == 0){
      ERr <- NULL
      E <- 0                                                    #mz
    } else {
      ERr <- t(sapply(ERr, constraint_to_row, varnames = varnames))
      E <- nrow(ERr)                                            #mz
      #colnames(ERr)[ncol(ERr)]<- "="                           #mz
    }
    if(length(IRr) == 0){
      IRr <- NULL
    } else {
      IRr <- t(sapply(IRr, constraint_to_row, varnames = varnames))
      #colnames(IRr)[ncol(IRr)]<- ">"                           #mz
    }

    if(E==0){
      R <- IRr[,-ncol(IRr)]}else{                               #mz
        if(length(IRr) == 0){R <- ERr[,-ncol(ERr)]}else{        #mz
          R <- rbind(ERr[,-ncol(ERr)],IRr[,-ncol(IRr)])}}       #mz
    if(is.vector(R)==TRUE){R<-matrix(R,nrow=1)}
    #mz
    r <- c(ERr[,ncol(ERr)],IRr[,ncol(IRr)])                     #mz
    list(Amat = R, difmin = r, E=E)                             #mz

  }), recursive = FALSE)

  #names(hyp_list) <- paste0(names(hyp_list), rep(1:(length(hyp_list)/2), each = 2))     #mz
  hyp_list
}






#suport functions matrices ####
parse_hypothesis <- function(varnames, hyp){
  names_est <- varnames
  # Clean varnames and hyp, to turn parameter names into legal object names
  # Might be better to do this elsewhere
  hyp <- rename_function(hyp)
  varnames <- rename_function(varnames)

  # Check if varnames occur in hyp.
  hyp_params <- params_in_hyp(hyp)

  # Find full varnames using partial matching
  match_names <- charmatch(hyp_params, varnames)
  # Any varnames that do not match?
  if(anyNA(match_names)){
    stop("Some of the parameters referred to in the 'hypothesis' do not correspond to parameter names of object 'x'.\n  The following parameter names in the 'hypothesis' did not match any parameters in 'x': ",
         paste(hyp_params[is.na(match_names)], collapse = ", "),
         "\n  The parameters in object 'x' are named: ",
         paste(varnames, collapse = ", "))
  }
  if(any(match_names == 0)){
    stop("Some of the parameters referred to in the 'hypothesis' matched multiple parameter names of object 'x'.\n  The following parameter names in the 'hypothesis' matched multiple parameters in 'x': ",
         paste(hyp_params[match_names == 0], collapse = ", "),
         "\n  The parameters in object 'x' are named: ",
         paste(varnames, collapse = ", "))
  }
  # Substitute partial names in hypothesis with full names
  for(par in 1:length(hyp_params)){
    hyp <- gsub(hyp_params[par], varnames[match_names[par]], hyp)
  }
  # Substitute parameter names with full names
  hyp_params <- varnames[match_names]

  # Start parsing hypothesis here

  legal_varnames <- sapply(hyp_params, grepl, pattern = "^[a-zA-Z\\.][a-zA-Z0-9\\._]{0,}$")
  if(!all(legal_varnames)){
    stop("Could not parse the names of the 'estimates' supplied to bain(). Estimate names must start with a letter or period (.), and can be a combination of letters, digits, period and underscore (_).\nThe estimates violating these rules were originally named: ", paste("'", names_est[!legal_varnames], "'", sep = "", collapse = ", "), ".\nAfter parsing by bain, these parameters are named: ", paste("'", hyp_params[!legal_varnames], "'", sep = "", collapse = ", "), call. = FALSE)
  }

  # Currently disabled: Is it a problem if there are parameters that don't occur in the hypothesis?
  # if(FALSE & any(!varnames %in% hyp_params)){
  #   message(
  #     "Some of the parameters in your model are not referred to in the 'hypothesis'. Make sure that this is what you intended.\n  Your hypothesis is: ",
  #     hyp,
  #     "\n  The parameters that are not mentioned are: ",
  #     paste(varnames[which(!varnames %in% hyp_params)], collapse = ", ")
  #   )
  # }
  # Ultimately, it would be best to check the reverse as well. This would also
  # allow us to make a string like "(a|b|c)" with the varnames, and insert it
  # into the regular expressions below. That would be faster than the complex
  # expressions currently used, and less prone to breakage.
  # varnames_in_hyp <- sapply(varnames, grepl, x = hyp)
  # if (FALSE) {
  #
  # }
  # End check
  hyp <- gsub("\\s", "", hyp)
  if(grepl("[><=]{2,}", hyp)) stop("Do not use combined comparison signs e.g., '>=' or '=='")

  original_hypothesis <- hyp_list <- strsplit(hyp, ";")[[1]]
  hyp_list <- lapply(hyp_list, function(x){ strsplit(x, "&")[[1]]})
  hyp_list <- lapply(hyp_list, function(x){unlist(lapply(x, expand_compound_constraints))})
  hyp_list <- lapply(hyp_list, function(x){unlist(lapply(x, expand_parentheses))})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, flip_inequality)})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, constraint_to_equation)})
  hyp_list <- lapply(hyp_list, function(x){sapply(x, order_terms)})

  n_constraints <- as.vector(sapply(hyp_list, function(x){c(sum(grepl("=", x)), sum(grepl("[<>]", x)))}))

  hyp_mat <- do.call(rbind, lapply(1:length(hyp_list), function(i){
    if(n_constraints[((i-1)*2)+1] == 0){
      ERr <- NULL
    } else {
      ERr <- t(sapply(hyp_list[[i]][grep("=", hyp_list[[i]])],
                      constraint_to_row, varnames = varnames))
    }
    if(n_constraints[((i-1)*2)+2] == 0){
      IRr <- NULL
    } else {
      IRr <- t(sapply(hyp_list[[i]][grep("[<>]", hyp_list[[i]])],
                      constraint_to_row, varnames = varnames))
    }
    rbind(ERr, IRr)
  }))

  list(hyp_mat = hyp_mat, n_constraints = n_constraints, original_hypothesis = original_hypothesis)
}

#' Expand compound constraints
#'
#' Takes a compound BAIN constraint, with multiple (in)equality operators, and
#' expands it into simple constraints.
#'
#'
#' @param hyp Character. A BAIN (in)equality constraint
#' @return A character vector with one element for each simple constraint
#' @keywords internal
expand_compound_constraints <- function(hyp){
  equality_operators <- gregexpr("[=<>]", hyp)[[1]]
  if(length(equality_operators) > 1){
    string_positions <- c(0, equality_operators, nchar(hyp)+1)
    return(sapply(1:(length(string_positions)-2), function(pos){
      substring(hyp, (string_positions[pos]+1), (string_positions[pos+2]-1))
    }))
  } else {
    return(hyp)
  }
}





#' Expand parentheses
#'
#' Takes a BAIN constraint with parentheses containing a "vector" of parameter
#' labels, and recursively expands the parenthesized "vectors".
#'
#'
#' @param hyp Character. A BAIN (in)equality constraint
#' @return A character vector with one element for each simple constraint
#' @keywords internal
expand_parentheses <- function(hyp){
  parenth_locations <- gregexpr("[\\(\\)]", hyp)[[1]]
  if(!parenth_locations[1] == -1){
    if(length(parenth_locations) %% 2 > 0) stop("Not all opening parentheses are matched by a closing parenthesis, or vice versa.")
    expanded_contents <- strsplit(substring(hyp, (parenth_locations[1]+1), (parenth_locations[2]-1)), ",")[[1]]
    if(length(parenth_locations) == 2){
      return(paste0(substring(hyp, 1, (parenth_locations[1]-1)), expanded_contents, substring(hyp, (parenth_locations[2]+1), nchar(hyp))))
    } else {
      return(apply(
        expand.grid(expanded_contents, expand_parentheses(substring(hyp, (parenth_locations[2]+1), nchar(hyp)))),
        1, paste, collapse = ""))
    }
  } else {
    return(hyp)
  }
}





#' Flip inequality
#'
#' Takes a BAIN constraint and flips elements on both sides of any "<"
#' operators so that only the ">" operator is used to define inequality
#' constraints.
#'
#'
#' @param hyp Character. A BAIN (in)equality constraint
#' @return Character
#' @keywords internal
flip_inequality <- function(hyp){
  if(grepl("<", hyp)){
    loc <- gregexpr("<", hyp)[[1]][1]
    return(paste0(substring(hyp, (loc+1)), ">", substring(hyp, 1, (loc-1))))
  } else {
    return(hyp)
  }
}





#' Constraint to equation
#'
#' Formats a BAIN constraint as an equation that can be evaluated. Adds scalars
#' to all parameters in the constraint, and adds a XXXconstant parameter to any
#' constants. Also adds "*" and "+" operators where necessary for evaluation.
#'
#'
#' @param hyp Character. A BAIN (in)equality constraint
#' @return Character
#' @keywords internal
constraint_to_equation <- function(hyp){
  # If scalar comes after variable name, move it in front
  hyp <- gsub("([\\.a-zA-Z][a-zA-Z0-9_\\.]{0,})\\*(\\d+)", "\\2*\\1", hyp, perl = TRUE)

  # When the string starts with a word, OR when a word is not preceded by
  # a number or *-sign, replace the "word" with "1*word"
  # hyp <- gsub("(^|(?<![\\*\\d]))([a-zA-Z][a-zA-Z0-9_]{0,})", "1*\\2", hyp, perl = TRUE)
  # Gu: add if not being preceded by any of .a-zA-Z0-9 , then ignore. otherwise x.A becomes x.1*A
  hyp <- gsub("(^|(?<![\\*\\d\\.a-zA-Z0-9_]))([a-zA-Z\\.][a-zA-Z0-9_\\.]{0,})", "1*\\2", hyp, perl = TRUE)

  # If a number starts with a period, add a leading zero
  # hyp <- gsub("(^|(?<!\\d))(\\.[0-9]+)", "0\\2", hyp, perl = TRUE)
  # Gu: add if being preceded by any of a-zA-Z0-9_, then ignore. Otherwise x.1 becomes x0.1
  hyp <- gsub("(^|(?<![a-zA-Z0-9_]))(\\.[0-9]+)", "0\\2", hyp, perl = TRUE)

  # When a number is directly followed by a word, insert a * sign
  # hyp <- gsub("(\\d)(?=[a-zA-Z][a-zA-Z0-9_]{0,})", "\\1*", hyp, perl = TRUE)
  # Gu: but not precededed by any of a-zA-Z_. Otherwise, x1x becomes x1*x
  hyp <- gsub("(?<![a-zA-Z_\\.])(\\d)(?=[a-zA-Z][a-zA-Z0-9_]{0,})", "\\1*", hyp, perl = TRUE) #####

  # When the string starts with a floating point number, OR
  # when a floating point number is not preceded by + or -, add a +
  # hyp <- gsub("(^|(?<![+-]))([0-9]{1,}\\.[0-9]{0,})", "+\\2", hyp, perl = TRUE)
  # Gu: add a condition: followed by a *, because otherwise x0.1x becomes x+0.1x.
  hyp <- gsub("(^|(?<![+-]))([0-9]{1,}\\.[0-9]{0,}(?=\\*))", "+\\2", hyp, perl = TRUE)

  # When the string starts with an integer number, OR
  # when an integer number is not preceded by + or -, add a +
  # hyp <- gsub("(^|(?<![\\.0-9+-]))([0-9]{1,})(?!\\.)", "+\\2", hyp, perl = TRUE)
  # Gu: add a condition: followed by a *, because otherwise x1x becomes x+1x.
  hyp <- gsub("(^|(?<![\\.0-9+-]))([0-9]{1,})(?!\\.)(?=\\*)", "+\\2", hyp, perl = TRUE)

  # Gu: To complement the above two, for a single number (constant), add a +
  hyp <- gsub("(?<=[=<>]|^)(\\d+|(\\d+\\.\\d+))(?=[=+-<>]|$)", "+\\1", hyp, perl = TRUE)

  # When a number is followed by =, >, <, +, or -, or is the last character of
  # the string, add "*XXXconstant" to the number
  # hyp <- gsub("(\\d)((?=[=<>+-])|$)", "\\1*XXXconstant", hyp, perl = TRUE)
  # Gu: but not preceded by word, number, period and  underline
  gsub("(?<![a-zA-Z0-9_\\.])(\\d+|(\\d+\\.\\d+))((?=[=<>+-])|$)", "\\1*XXXconstant", hyp, perl = TRUE)
}






#' Order terms
#'
#' Moves all terms on the right hand side of a BAIN constraint, which has been
#' formatted as an equation, to the left hand side.
#'
#'
#' @param hyp Character. A BAIN (in)equality constraint formatted as equation.
#' @return Character
#' @keywords internal
order_terms <- function(hyp){
  eq_location <- gregexpr("[=<>]", hyp)[[1]]
  rhs <- substring(hyp, eq_location+1, nchar(hyp))
  rhs <- gsub("(?=[+-])", "-", rhs, perl = TRUE)
  paste0(rhs, substring(hyp, 1, eq_location))
}





#' Constraint to row
#'
#' Evaluate a BAIN constraint, which has been formatted as an equation, with
#' all terms moved to the left hand side, and return a single row of a BAIN
#' (in)equality constraint matrix.
#'
#'
#' @param hyp Character. A BAIN (in)equality constraint formatted as equation,
#' with all terms on the left hand side.
#' @return Numeric vector.
#' @keywords internal
constraint_to_row <- function(varnames, hyp){
  e <- new.env()
  objects <- c(varnames, "XXXconstant")
  invisible(sapply(objects, assign, value = 0, envir = e))
  constraint_expression <- parse(text = substring(hyp, 1, nchar(hyp)-1))
  e[["XXXconstant"]] <- -1
  equal_to <- eval(constraint_expression, envir = e)
  e[["XXXconstant"]] <- 0

  c(sapply(varnames, function(x){
    e[[x]] <- 1
    constant <- eval(constraint_expression, envir = e)
    e[[x]] <- 0
    constant
  }), equal_to)
}

params_in_hyp <- function(hyp){
  params_in_hyp <- trimws(unique(strsplit(hyp, split = "[ =<>,\\(\\);&\\*+-]+", perl = TRUE)[[1]]))
  params_in_hyp <- params_in_hyp[!sapply(params_in_hyp, grepl, pattern = "^[0-9]*\\.?[0-9]+$")]
  params_in_hyp[grepl("^[a-zA-Z]", params_in_hyp)]
}

#' @importFrom utils tail
rename_function <- function(text){
  fulltext <- paste(text, collapse = "")
  new_names <- names_est <- text
  #if(grepl("[\\(\\)]", fulltext)){
  #  text <- gsub("\\(", "___O___", text)
  #  text <- gsub("\\)", "___C___", text)
  #}
  text[text == "(Intercept)"] <- "Intercept"
  if(grepl(":", fulltext)){
    text <- gsub(":", "___X___", text)
  }

  if(grepl("mean of ", fulltext)){
    text <- gsub("mean of the differences", "difference", text)
    text <- gsub("mean of ", "", text)
  }

  # If any variables are subsetted from data.frames: remode the df part of the name
  remove_df <- sapply(text, grepl, pattern = "[\\]\\$]+", perl = TRUE)
  if(any(remove_df)){
    text[remove_df] <- sapply(text[remove_df], function(x){
      tmp_split <- strsplit(x, "[\\]\\$]+", perl = TRUE)[[1]]
      if(length(tmp_split)==1){
        x
      } else {
        tail(tmp_split, 1)
      }
    })
  }

  text
}

#' @importFrom utils tail
rename_estimate <- function(estimate){

  new_names <- names_est <- names(estimate)
  if(any(new_names == "(Intercept)")) new_names[match(new_names, "(Intercept)")] <- "Intercept"
  if(is.null(names_est)){
    stop("The 'estimates' supplied to bain() were unnamed. This is not allowed, because estimates are referred to by name in the 'hypothesis' argument. Please name your estimates.")
  }

  if(length(new_names) < 3){
    new_names <- gsub("mean of the differences", "difference", new_names)
    new_names <- gsub("mean of ", "", new_names)
  }

  # If any variables are subsetted from data.frames: remode the df part of the name
  remove_df <- sapply(new_names, grepl, pattern = "[\\]\\$]+", perl = TRUE)
  if(any(remove_df)){
    new_names[remove_df] <- sapply(new_names[remove_df], function(x){
      tmp_split <- strsplit(x, "[\\]\\$]+", perl = TRUE)[[1]]
      if(length(tmp_split)==1){
        x
      } else {
        tail(tmp_split, 1)
      }
    })
  }

  # Any interaction terms: replace : with _X_
  new_names <- gsub(":", "___X___", new_names)

  legal_varnames <- sapply(new_names, grepl, pattern = "^[a-zA-Z\\.][a-zA-Z0-9\\._]{0,}$")
  if(!all(legal_varnames)){
    stop("Could not parse the names of the 'estimates' supplied to bain(). Estimate names must start with a letter or period (.), and can be a combination of letters, digits, period and underscore (_).\nThe estimates violating these rules were originally named: ", paste("'", names_est[!legal_varnames], "'", sep = "", collapse = ", "), ".\nAfter parsing by bain, these parameters are named: ", paste("'", new_names[!legal_varnames], "'", sep = "", collapse = ", "), call. = FALSE)
  }
  names(estimate) <- new_names
  estimate
}


