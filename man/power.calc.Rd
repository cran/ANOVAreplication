\name{power.calc}
\alias{power.calc}
\title{
Power Calculator}
\description{
Calculates the power for the prior predictive check against equal means}
\usage{
power.calc(n.r,posterior,g.m,p.sd,statistic,Amat=0L,exact=0L,difmin=0L,effectsize=FALSE,
           alpha=.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.r}{vector with the sample size per group (i.e., n_{jr}) for new study (i.e., y_r).}
  \item{posterior}{matrix (e.g., the output of Gibbs.ANOVA) with samples from the posterior based on the original data (i.e., y_o).}
  \item{g.m}{integer; the population value for the equal means in the alternative distribution. We advice to specify the grand mean of the study variables in the original dataset.}
  \item{p.sd}{integer; the population value for the pooled standard deviation in the alternative distribution. We advice to specify the pooled standard deviation for the study variables in the original dataset.}
  \item{statistic}{the type of hypothesis to be evaluated: "ineq" for inequality constrained means, "dif" for inequality constraints plus minimum differences between means, "exact" for specific values for the means.}
  \item{Amat}{p by q matrix, where p is the number of means in the ANOVA model, and q is the number of constraints to be imposed on the model. Each row represents one constraint where the parameter with the lower value according to the constraint receives the value -1, and the parameter with the higher value according to the constraint receives the value 1. Other parameters within the same row obtain the value 0.}
\item{exact}{vector of length p, where p is the number of means in the ANOVA model, with the exact values of the constrained hypothesis.}
\item{difmin}{vector of length q with the minimum difference per constraint as specified in \code{Amat}.}
\item{effectsize}{logical; If TRUE the values in \code{difmin} are interpreted as Cohen's d.}
\item{alpha}{integer; the level of alpha that should be taken into account while calculating the required sample size.}

}
\value{
\item{power}{The acquired power given the input}
\item{rejection.value}{The 1-alpha'th percentile of the null distribution. The proportion of H1 larger than this value constitues power.}
}
\references{
Zondervan-Zwijnenburg, M.A.J., Van de Schoot, R., & Hoijtink, H. (2017). Testing ANOVA replication by means of the prior predictive p-value. }
\author{
M. A. J. Zondervan-Zwijnenburg}
\seealso{
See also \code{\link{runShiny}}, \code{\link{Gibbs.ANOVA}}, \code{\link{Fbar.ineq}}, \code{\link{Fbar.dif}}, and \code{\link{Fbar.exact}}, \code{\link{prior.predictive.check}}, \code{\link{sample.size.calc}}.
}
\examples{
#analysis original data
data_o <- data.frame(y=ChickWeight$weight,g=ChickWeight$Diet)
g.m <- mean(data_o$y)
#compute pooled sd
sd.g <- aggregate(data_o$y,by=list(data_o$g),sd)[,2]
n.g <- table(data_o$g)
p.sd <- ((n.g[1]-1)*sd.g[1]+(n.g[2]-1)*sd.g[2]+(n.g[3]-1)*sd.g[3]+(n.g[4]-1)*sd.g[4])/(sum(n.g)-4)

means <- aggregate(data_o$y,by=list(data_o$g),mean)[,2]

Gibbs.ANOVA(data_o,it=200,burnin=50) #we advise >1000 iterations, >500 burnin for reliable results

power.calc(n.r=c(20,21,22,23),posterior=output_m,g.m=g.m,p.sd=p.sd,
           statistic="exact",exact=means,alpha=.05)

}
\keyword{htest}
