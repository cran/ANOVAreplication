library(shiny)
source("helpers.R")

shinyServer(function(input, output, session) {

  #input original data ####
  data_o <- reactive({
    # input$file_o will be NULL initially. After the user selects and uploads a file, it will be a data frame with 'name', 'size', 'type', and 'datapath'
    # columns. The 'datapath' column will contain the local filenames where the data can be found.
    inFile <- input$file_o
    if (is.null(inFile)){return()}
    df_o <- read.csv(inFile$datapath, header=input$header, sep=input$sep)
    names(df_o) <- c("y","g")
    df_o
  })

  #input original data descriptives ####
  data_descr <- reactive({
    if (input$runButton_descr == 0){return()}
    p <- input$n1
    means_1 <- lapply(1:p, function(i) input[[paste0('groupmean', i)]])
    means_2 <- lapply(1:p, function(i) as.numeric(means_1))
    means <- do.call(c,means_2)

    sds_1 <- lapply(1:p, function(i) input[[paste0('groupsd', i)]])
    sds_2 <- lapply(1:p, function(i) as.numeric(sds_1))
    sds <- do.call(c,sds_2)

    n_1 <- lapply(1:p, function(i) input[[paste0('groupn', i)]])
    n_2 <- lapply(1:p, function(i) as.numeric(n_1))
    n <- do.call(c,n_2)

    y <- list(NA); group <- list(NA)
    for (i in 1:p){
      y[[i]] <- rnorm2(n[i],means[i],sds[i])
      group[[i]] <- rep(i,n[i])
    }
    y <- unlist(y)
    group <- unlist(group)
    #save the data for gibbs
    if (input$typepriorinput==2){
      data_o.gibbs <<- list(y=y,g=group,p=p,N=sum(n))}
    data.frame(y,group)
  })

  #tables
  output$contents_o <- renderTable({
    if (is.null(data_o())){return()}
    data_o.work <- data_o()
    data_o.work
  })
  output$contents_descr <- renderTable({
    if (is.null(data_descr())){return()}
    data_o.work <- data_descr()
    data_o.work
  })

  #summary stats own data
  #summary and data preparing for own data
  output$summary_o1a <- renderText({
    if (is.null(data_o())){return()}
    data_o.work <- data_o()

    data_o.p <- length(unique(data_o.work[,2])) #number of groups
    data_o.n <- length(data_o.work[,1])         #sample size
    #save the data for gibbs
    if (input$typepriorinput==1){
      data_o.gibbs <<- list(y=data_o.work[,1],g=data_o.work[,2],p=data_o.p,N=data_o.n)}

    paste("The total sample size is:", data_o.n)
  })
  output$summary_o1b <- renderPrint({
    if (is.null(data_o())){return(invisible())}
    table(data_o()$g) })
  output$summary_o1c <- renderPrint({
    if (is.null(data_o())){return(invisible())}
    setNames(aggregate(data_o()$y,by=list(data_o()$g),mean), c("Group","Mean"))})

  output$summary_o2a <- renderText({
    if (is.null(data_descr())){return()}
    data_o.work <- data_descr()

    data_o.p <- length(unique(data_o.work[,2])) #number of groups
    data_o.n <- length(data_o.work[,1])         #sample size
    #save the data for gibbs
    if (input$typepriorinput==2){
      data_o.gibbs <<- list(y=data_o.work[,1],g=data_o.work[,2],p=data_o.p,N=data_o.n)}

    paste("The total sample size is:", data_o.n)
  })
  output$summary_o2b <- renderPrint({
    if (is.null(data_descr())){return(invisible())}
    table(data_descr()$g) })
  output$summary_o2c <- renderPrint({
    if (is.null(data_descr())){return(invisible())}
    setNames(aggregate(data_descr()$y,by=list(data_descr()$g),mean), c("Group","Mean"))})

  #Bayesian analysis ####
  output$helptext1 <- renderUI({
    #if (input$runButton == 0){return()}
    if (is.null(output_Gibbs())){return(invisible())}

    withMathJax(
      helpText(
        "The table below gives the output of the Bayesian analysis on the original data.",
        "In the first column the means of the marginal posterior distributions are given, while the second column provides the standard deviation of these distributions.",
        "The third and fourth column show the limits of the 95% credible interval."))
  })

  output_Gibbs <- reactive({
    if (input$runButton == 0){return()}

    isolate({
      #R program for Gibbs sampling from full conditionals in OLS example
      #number of iterations
      seed_G <- as.numeric(input$seed_G)
      if(seed_G!=0){
        set.seed(seed=seed_G)}
      burnin=500
      I<-as.numeric(input$n.iter)+as.numeric(input$n.burnin)
      #read only observations with complete information
      data<- data_o.gibbs
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
      par1<<-cbind(b1,s1)[-c(1:burnin),] #100 for burnin
      par2<<-cbind(b2,s2)[-c(1:burnin),]
      output_G <- rbind(par1,par2)
      colnames(output_G) <- paste("Mean",1:(G+1))
      colnames(output_G)[G+1] <- "SD"
      #output_G <<- output_G
      output_G
    })
  })


  output$posterior<- renderTable({
    if (is.null(output_Gibbs())){return(invisible())}

    output_G <- output_Gibbs()

    posterior <- cbind(round(apply(output_G,2,mean),4),
                       round(apply(output_G,2,sd),4),
                       round(apply(output_G,2,quantile,.05),3),
                       round(apply(output_G,2,quantile,.95),3))
    colnames(posterior) <- c("Estimate","SD","5%","95%")

    p <- data_o.gibbs$p
    rownames(posterior) <- paste("Mean",1:(p+1))
    rownames(posterior)[(p+1)] <- "SD"
    posterior
  })

  output$helptext0 <- renderText({
    if (is.null(output_Gibbs())){return(invisible())}
    #if (input$runButton == 0){return()}
    paste(
      "The plots below help you to examine the convergence of your model. The right column shows traceplots.",
      "The traceplots should look like fat catterpillars, that is, they should have stable means and variances.",
      "If the traceplots do not have stable means and variances, you should increase the number of iterations to obtain the posterior distribution.",
      "The plots on the left are the posterior distributions: the outcomes of the analysis.")
  })

  output$convergence<- renderPlot({
    if (is.null(output_Gibbs())){return(invisible())}
    output_G <- output_Gibbs()

    par(mfrow = c(data_o.gibbs$p+1, 2), mar = rep(2, 4))

    for (j in 1:(data_o.gibbs$p+1)) {

      hist(output_G[,j], breaks="Scott", main=colnames(output_G)[j], xlab="Value")
      abline(v=quantile(output_G[,j], probs=c(0.025, 0.975)), col="blue")

      plot.ts(cbind(par1[,j],par2[,j]),col=c("blue","green"), plot.type="single",ylab="")

    }

  })

  data_r <- reactive({
    inFile1 <- input$file_r
    if (is.null(inFile1)){return()}
    df_r <- read.csv(inFile1$datapath, header=input$header_r, sep=input$sep_r)
    names(df_r) <- c("y","g")
    df_r
  })


  #input replication data descriptives ####
  data_descr_r <- reactive({
    if (input$runButton_descr_r == 0){return()}
    p <- input$nr1
    means_1 <- lapply(1:p, function(i) input[[paste0('groupmean_r', i)]])
    means_2 <- lapply(1:p, function(i) as.numeric(means_1))
    means <- do.call(c,means_2)

    sds_1 <- lapply(1:p, function(i) input[[paste0('groupsd_r', i)]])
    sds_2 <- lapply(1:p, function(i) as.numeric(sds_1))
    sds <- do.call(c,sds_2)

    n_1 <- lapply(1:p, function(i) input[[paste0('groupn_r', i)]])
    n_2 <- lapply(1:p, function(i) as.numeric(n_1))
    n <- do.call(c,n_2)

    y <- list(NA); group <- list(NA)
    for (i in 1:p){
      y[[i]] <- rnorm2(n[i],means[i],sds[i])
      group[[i]] <- rep(i,n[i])
    }
    y <- unlist(y)
    group <- unlist(group)
    data.frame(y,group)
  })

  #tables
  output$contents_r <- renderTable({
    if (is.null(data_r())){return()}
    data_r.work <<- data_r()
    data_r.work
  })
  output$contents_descr_r <- renderTable({
    if (is.null(data_descr_r())){return()}
    data_r.work <<- data_descr_r()
    data_r.work
  })

  observe({
    dataRpresent <- !is.null(data_descr_r()) | !is.null(data_r())
    # Change the selected tab.
    if (dataRpresent) {
      updateTabsetPanel(session, "outputTabset", selected = "Replication Data")}
  })

  output$summary_r1a <- renderText({
    if (is.null(data_r())){return()}
    data_r.work <- data_r()
    data_r.n <- length(data_r.work[,1])         #sample size
    paste("The total sample size is:", data_r.n)
  })
  output$summary_r1b <- renderPrint({
    if (is.null(data_r())){return(invisible())}
    table(data_r()$g) })
  output$summary_r1c <- renderPrint({
    if (is.null(data_r())){return(invisible())}
    setNames(aggregate(data_r()$y,by=list(data_r()$g),mean), c("Group","Mean"))})

  #summary stats generated rep data
  output$summary_r2a <- renderText({
    if (is.null(data_descr_r())){return()}
    data_r.work <- data_descr_r()

    data_r.p <- length(unique(data_r.work[,2])) #number of groups
    data_r.n <- length(data_r.work[,1])         #sample size
    #save the data for gibbs

    paste("The total sample size is:", data_r.n)
  })
  output$summary_r2b <- renderPrint({
    if (is.null(data_descr_r())){return(invisible())}
    table(data_descr_r()$g) })
  output$summary_r2c <- renderPrint({
    if (is.null(data_descr_r())){return(invisible())}
    setNames(aggregate(data_descr_r()$y,by=list(data_descr_r()$g),mean), c("Group","Mean"))})

  observe({
    # store means of original data (after posterior is observed)
    if (is.null(output_Gibbs())){return(invisible())}
    x <- round(aggregate(data_o.gibbs$y,by=list(data_o.gibbs$g),mean)[,2],3)

    # use the means as prefilled input for exact values in ppc and sampcalc
    updateTextInput(session, "exactval", value = paste(x))
    updateTextInput(session, "exactval_N", value = paste(x))
  })



  #prior predictive check ####
  results.ppc <- reactive({
    if (input$runButton_priorpred == 0){return()}
    validate(need(output_Gibbs(),message="Obtain the posterior distribution in the Original Study tab first."))
    validate(need(!is.null(data_r())|!is.null(data_descr_r()),
                  message="Upload data or provide descriptives in the New Study tab first."))

    isolate({
      if (input$typepriorinput == 1){data_o.work <- data_o()}
      if (input$typepriorinput == 2){data_o.work <- data_descr()}
      if (input$typerepinput == 1){data_r.work <- data_r()}
      if (input$typerepinput == 2){data_r.work <- data_descr_r()}

      n.r <- as.numeric(table(data_r.work[,2]))
      if (input$typehypothesis == "ineq"){

        Amat_1a <- lapply(1:input$nh, function(i) input[[paste0('Amat_a', i)]])
        Amat_2a <- lapply(1:input$nh,function(i) unlist(Amat_1a[[i]]))
        Amat_3a <- do.call(c,Amat_2a)

        Amat_1b <- lapply(1:input$nh, function(i) input[[paste0('Amat_b', i)]])
        Amat_2b <- lapply(1:input$nh,function(i) unlist(Amat_1b[[i]]))
        Amat_3b <- do.call(c,Amat_2b)

        Amat <- matrix(0,nrow = input$nh,ncol=data_o.gibbs$p)
        for(i in 1:input$nh){
          Amat[i,Amat_3a[i]] <- 1
          Amat[i,Amat_3b[i]] <- -1}
        validate(need(Amat,message="Provide inequality constraints."))
        validate(need(all(apply(Amat,1,sum)==0),message="Parameter cannot be > itself, adjust constraints"))

        exact=0L

        if (input$addDif==0){difmin=0L
        Fstat <<- Fbar.ineq(data_r.work,Amat)[1]
        statistic="ineq"; effectsize=FALSE}
        if (input$addDif==1){ #absolute differences

          difmin_1 <- lapply(1:input$nh, function(i) input[[paste0('difmin', i)]])
          difmin_2 <- lapply(1:input$nh,function(i) unlist(strsplit(difmin_1[[i]],",")))
          difmin_3 <- lapply(1:input$nh,function(i) as.numeric(difmin_2[[i]]))
          difmin <- do.call(c,difmin_3)

          validate(need(difmin!="",message="Provide minimum differences."))

          Fstat <<- Fbar.dif(data_r.work,Amat,difmin)[1]
          statistic="dif"; effectsize=FALSE}
        if (input$addDif==2){ #effect size differences

          difmin_1 <- lapply(1:input$nh, function(i) input[[paste0('difmin', i)]])
          difmin_2 <- lapply(1:input$nh,function(i) unlist(strsplit(difmin_1[[i]],",")))
          difmin_3 <- lapply(1:input$nh,function(i) as.numeric(difmin_2[[i]]))
          difmin <- do.call(c,difmin_3)

          validate(need(difmin!="",message="Provide minimum differences."))

          Fstat <<- Fbar.dif(data_r.work,Amat,difmin,effectsize=TRUE)[1]
          statistic="dif"; effectsize=TRUE}
      }

      if (input$typehypothesis == "exact"){
        Amat = 0L; difmin=0L
        exact = as.numeric(unlist(strsplit(input$exactval,",")))
        validate(need(exact!="",message="Provide the values to replicate."))
        Fstat <<- Fbar.exact(data_r.work,exact)[1]
        statistic="exact"; effectsize=FALSE}

      output_G <- output_Gibbs()
      outputppc <<- prior.predictive.check(n=n.r,posterior=output_G,F_obs=Fstat,statistic=statistic,
                                           Amat=Amat,exact=exact,difmin=difmin,effectsize=effectsize,seed=as.numeric(input$seed))
      Fbar <<- Fps
    })
  })

  output$priorpred <- renderPrint({
    if (is.null(results.ppc())){return(invisible())}
    outputppc
  })

  #ppc plot ####
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]}

  observe({
    if (is.null(results.ppc())){return()}
    #par(mar = rep(2, 4))
    h <<- hist(Fbar[-which(Fbar==Mode(Fbar))],breaks=20)
  })

  plotInput <- function(){
    if(sum(Fbar==Mode(Fbar))/length(Fbar)>=.10){
      hist(Fbar[-which(Fbar==Mode(Fbar))],ylim=c(0,max(sum(Fbar==Mode(Fbar)),h$counts[1])),
           xlim=c(0,max(max(Fbar),Fstat)),
           xlab=expression(italic(bar(F)[bold(y)])),ylab="Frequency",main="",breaks=20)
      segments(x0=Mode(Fbar),y0=0,x1=Mode(Fbar),y1=sum(Fbar==Mode(Fbar)),col="black",lwd=5)
      abline(v=Fstat,col="red")
    }else{
      hist(Fbar,freq=TRUE,breaks=seq(0,max(Fbar),length.out=40),xlim=c(0,max(max(Fbar),Fstat)),
           xlab=expression(italic(bar(F)[bold(y)])),ylab="Frequency",main="")
      abline(v=Fstat,col="red")}
  }

  output$downloadPlot <- downloadHandler(
    filename="ANOVAreplicationTest.pdf",
    content <- function(file) {
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      pdf(file,width=6,height=4)
      plotInput()
      dev.off()

      tempPlot <- file.path(tempdir(), "ANOVAreplicationTest.pdf")
      file.copy("ANOVAreplicationTest.pdf", tempPlot, overwrite = TRUE)
    }
  )

  #plot ####
  output$Fps<- renderPlot({
    if (is.null(results.ppc())){return()}

    Mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]}

    if(sum(Fbar==Mode(Fbar))/length(Fbar)>=.10){
      h <- hist(Fbar[-which(Fbar==Mode(Fbar))],ylim=c(0,sum(Fbar==Mode(Fbar))),xlab=expression(italic(bar(F)[bold(y)][r])),ylab="Frequency",main="",breaks=20)
      hist(Fbar[-which(Fbar==Mode(Fbar))],ylim=c(0,max(sum(Fbar==Mode(Fbar)),h$counts[1])),
           xlim=c(0,max(max(Fbar),Fstat)),
           xlab=expression(italic(bar(F)[bold(y)])),ylab="Frequency",main="",breaks=20)
      segments(x0=Mode(Fbar),y0=0,x1=Mode(Fbar),y1=sum(Fbar==Mode(Fbar)),col="black",lwd=5)
      abline(v=Fstat,col="red")
    }else{
      hist(Fbar,freq=TRUE,breaks=seq(0,max(Fbar),length.out=40),xlim=c(0,max(max(Fbar),Fstat)),
           xlab=expression(italic(bar(F)[bold(y)])),ylab="Frequency",main="")
      abline(v=Fstat,col="red")}
  })

  output$helptext2 <- renderText({
    if (is.null(results.ppc())){return()}
    paste(
      "Below you can find the results of the prior predictive check. The summary of the distribution of F-bar given the original data",
      " (also plotted in the histogram) demonstrates likely F-bar values for future studies given the parameter estimates from the original study.",
      "'F-bar replication data' is the F-bar value obtained in the replication study (indicated with a red line in the histogram). ",
      "The prior predictive p-value indicates how extreme the obtained F-bar in the replication study is compared to the F-bar values that are expected based on the original study.")
  })
  output$downloadData <- downloadHandler(filename="ANOVAreplicationTest.csv",
                                         content = function(file){write.csv(Fbar,file,row.names=FALSE)}
  )

  observe({
    # Change the selected tab.
    if (input$runButton_priorpred == 1) {
      updateTabsetPanel(session, "outputTabset", selected = "Replication Test Results")}
  })


  #sample size calculator ####
  results.ppp <- reactive({
    if (input$runButton_sampcalc == 0){return()}
    validate(need(output_Gibbs(),message="Obtain the posterior distribution in the Original Study tab first."))

    isolate({
      if (input$typepriorinput == 1){data_o.work <- data_o()}
      if (input$typepriorinput == 2){data_o.work <- data_descr()}

      if (input$typehypothesis_N == "ineq"){
        Amat_1a <- lapply(1:input$nh_N, function(i) input[[paste0('Amat_N_a', i)]])
        Amat_2a <- lapply(1:input$nh_N,function(i) unlist(Amat_1a[[i]]))
        Amat_3a <- do.call(c,Amat_2a)

        Amat_1b <- lapply(1:input$nh_N, function(i) input[[paste0('Amat_N_b', i)]])
        Amat_2b <- lapply(1:input$nh_N,function(i) unlist(Amat_1b[[i]]))
        Amat_3b <- do.call(c,Amat_2b)

        Amat_N <- matrix(0,nrow = input$nh_N,ncol=data_o.gibbs$p)
        for(i in 1:input$nh_N){
          Amat_N[i,Amat_3a[i]] <- 1
          Amat_N[i,Amat_3b[i]] <- -1}

        exact_N=0L;

        if (input$addDif_N==0){difmin_N=0L; effectsize_N = FALSE; statistic="ineq"}
        if (input$addDif_N==1){ #absolute differences
          difmin_1 <- lapply(1:input$nh_N, function(i) input[[paste0('difmin_N', i)]])
          difmin_2 <- lapply(1:input$nh_N,function(i) unlist(strsplit(difmin_1[[i]],",")))
          difmin_3 <- lapply(1:input$nh_N,function(i) as.numeric(difmin_2[[i]]))
          difmin_N <- do.call(c,difmin_3); effectsize_N = FALSE; statistic="dif"}
        if (input$addDif_N==2){ #effect size differences
          difmin_1 <- lapply(1:input$nh_N, function(i) input[[paste0('difmin_N', i)]])
          difmin_2 <- lapply(1:input$nh_N,function(i) unlist(strsplit(difmin_1[[i]],",")))
          difmin_3 <- lapply(1:input$nh_N,function(i) as.numeric(difmin_2[[i]]))
          difmin_N <- do.call(c,difmin_3)
          effectsize_N = TRUE; statistic="dif"}
      }
      if (input$typehypothesis_N == "exact"){
        Amat_N = 0L; difmin_N=0L; statistic="exact"; effectsize_N = FALSE
        exact_N = as.numeric(unlist(strsplit(input$exactval_N,",")))}

      output_G <- output_Gibbs()
      #nsample= dim(output_G)[1]
      #sample.r <- sample(1:dim(output_G)[1],nsample)

      outputsampcalc <<- power.ppp(start_n=as.numeric(input$start_n),itmax=as.numeric(input$maxit),
                                   powtarget=as.numeric(input$Powtarget),powmargin=as.numeric(input$Powmargin),
                                   posterior=output_G,g.m=mean(data_o.work$y),
                                   statistic=statistic,
                                   Amat=Amat_N,exact=exact_N,difmin=difmin_N,effectsize=effectsize_N,
                                   nmax=as.numeric(input$maxN),alpha=as.numeric(input$alpha))

    })
  })

  observe({
    # Change the selected tab.
    if (input$runButton_sampcalc == 1) {
      updateTabsetPanel(session, "outputTabset", selected = "Sample Size & Power Output")}
  })

  output$sampcalc <- renderPrint({
    if (is.null(results.ppp())){return(invisible())}
    outputsampcalc
  })

  output$Fpspower<- renderPlot({
    if (is.null(results.ppp())){return()}
    hist(Fps.power.H0,freq=FALSE,col=rgb(1,0,0,1/4),border=rgb(1,0,0,1/2),main="",xlab=expression(italic(bar(F)[bold(y)]))) #null with true effect
    hist(Fps.power.H1,freq=FALSE,col=rgb(0,0,1,1/4),border=rgb(0,0,1,1/2),add=TRUE) #H1, means equal
    abline(v=rej.value,col=rgb(1,0,1,1/2),lwd=2)
  })

  output$helptext3 <- renderText({
    if (is.null(results.ppp())){return()}
    paste(
      "Below you can find the results of the sample size calculator for the prior predictive check. ",
      "First the output provides the reason to stop the iterative sample size calculations. ",
      "Subsequently, the matrix with the output is given. ",
      "The number in column n per group is the sample size per group, the value on the right is the associated power.",
      "The histogram is based on information for the last sample size and power calculated.",
      "The red distribution is composed of F-bars for datasets from a population in which replication holds (i.e., the null distribution). ",
      "The blue distribution shows F-bars from a population with equal means for which replication should be rejected (i.e., the alternative distribution).",
      "The vertical line indicates the critical value located at 1-alpha'th percentile of the null distribution. ",
      "The proportion of the alternative distribution at the right side of the critical value constitutes statical power. ")
  })
  
  #power calculator ####
  results.basicpower <- reactive({
    if (input$runButton_powercalc == 0){return()}
    validate(need(output_Gibbs(),message="Obtain the posterior distribution in the Original Study tab first."))
    
    isolate({
      if (input$typepriorinput == 1){data_o.work <- data_o()}
      if (input$typepriorinput == 2){data_o.work <- data_descr()}
      
      n.r <- as.numeric(unlist(strsplit(input$nF,",")))
      
      if (input$typehypothesis_N == "ineq"){
        Amat_1a <- lapply(1:input$nh_N, function(i) input[[paste0('Amat_N_a', i)]])
        Amat_2a <- lapply(1:input$nh_N,function(i) unlist(Amat_1a[[i]]))
        Amat_3a <- do.call(c,Amat_2a)
        
        Amat_1b <- lapply(1:input$nh_N, function(i) input[[paste0('Amat_N_b', i)]])
        Amat_2b <- lapply(1:input$nh_N,function(i) unlist(Amat_1b[[i]]))
        Amat_3b <- do.call(c,Amat_2b)
        
        Amat_N <- matrix(0,nrow = input$nh_N,ncol=data_o.gibbs$p)
        for(i in 1:input$nh_N){
          Amat_N[i,Amat_3a[i]] <- 1
          Amat_N[i,Amat_3b[i]] <- -1}
        
        exact_N=0L;
        
        if (input$addDif_N==0){difmin_N=0L; effectsize_N = FALSE; statistic="ineq"}
        if (input$addDif_N==1){ #absolute differences
          difmin_1 <- lapply(1:input$nh_N, function(i) input[[paste0('difmin_N', i)]])
          difmin_2 <- lapply(1:input$nh_N,function(i) unlist(strsplit(difmin_1[[i]],",")))
          difmin_3 <- lapply(1:input$nh_N,function(i) as.numeric(difmin_2[[i]]))
          difmin_N <- do.call(c,difmin_3); effectsize_N = FALSE; statistic="dif"}
        if (input$addDif_N==2){ #effect size differences
          difmin_1 <- lapply(1:input$nh_N, function(i) input[[paste0('difmin_N', i)]])
          difmin_2 <- lapply(1:input$nh_N,function(i) unlist(strsplit(difmin_1[[i]],",")))
          difmin_3 <- lapply(1:input$nh_N,function(i) as.numeric(difmin_2[[i]]))
          difmin_N <- do.call(c,difmin_3)
          effectsize_N = TRUE; statistic="dif"}
      }
      if (input$typehypothesis_N == "exact"){
        Amat_N = 0L; difmin_N=0L; statistic="exact"; effectsize_N = FALSE
        exact_N = as.numeric(unlist(strsplit(input$exactval_N,",")))}
      
      output_G <- output_Gibbs()
      #nsample= dim(output_G)[1]
      #sample.r <- sample(1:dim(output_G)[1],nsample)
      
      outputpowercalc <<- power.basic(nF=n.r,posterior=output_G,g.m=mean(data_o.work$y),
                                   statistic=statistic,Amat=Amat_N,exact=exact_N,difmin=difmin_N,effectsize=effectsize_N,
                                   alpha=as.numeric(input$alpha))
      
    })
  })
  
  observe({
    # Change the selected tab.
    if (input$runButton_powercalc == 1) {
      updateTabsetPanel(session, "outputTabset", selected = "Sample Size & Power Output")}
  })
  
  output$powercalc <- renderPrint({
    if (is.null(results.basicpower())){return(invisible())}
    outputpowercalc
  })
  
  output$Fpsbasicpower<- renderPlot({
    if (is.null(results.basicpower())){return()}
    hist(Fps.power.H0,freq=FALSE,col=rgb(1,0,0,1/4),border=rgb(1,0,0,1/2),main="",xlab=expression(italic(bar(F)[bold(y)]))) #null with true effect
    hist(Fps.power.H1,freq=FALSE,col=rgb(0,0,1,1/4),border=rgb(0,0,1,1/2),add=TRUE) #H1, means equal
    abline(v=rej.value,col=rgb(1,0,1,1/2),lwd=2)
  })
  
  output$helptext3power <- renderText({
    if (is.null(results.basicpower())){return()}
    paste(
      "Below you can find the results of the power calculator for the prior predictive check. ",
      "In print, the observed power and the 1-alpha'th value of the null distribution are provided. ",
      "The red distribution is composed of F-bars for datasets from a population in which replication holds (i.e., the null distribution). ",
      "The blue distribution shows F-bars from a population with equal means for which replication should be rejected (i.e., the alternative distribution).",
      "The vertical line indicates the critical value located at 1-alpha'th percentile of the null distribution. ",
      "The proportion of the alternative distribution at the right side of the critical value constitutes statical power. ")
  })

})
