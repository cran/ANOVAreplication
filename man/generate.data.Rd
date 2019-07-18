\name{generate.data}
\alias{generate.data}
\title{Generate Data}
\description{
Generates data of sample size n with the exact specified mean and standard deviation.}
\usage{
generate.data(n,mean,sd)
}
\arguments{
  \item{n}{integer; sample size.}
  \item{mean}{integer; mean of the data.}
  \item{sd}{integer; standard deviation of the data.}
}
\examples{
#simple example
data1 <- generate.data(n=10,mean=5,sd=1)

#multiple independent groups
means = c(0.36,-0.19,-0.18)
sds = c(1.08,0.53,0.81)
n = c(28,28,28) #N = 84
y <- list(NA)
for (i in 1:3){y[[i]] <- generate.data(n[i],means[i],sds[i])}
y <- unlist(y)
group <- c(rep(0,n[1]),rep(1,n[2]),rep(2,n[3]))
p <- length(unique(group))
data=data.frame(y=y,g=group)
aggregate(data$y,by=list(data$g),mean)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen }% use one of  RShowDoc("KEYWORDS")
