\name{prior.predictive.check}
\alias{prior.predictive.check}
\title{Prior predictive check}
\description{
Uses the prior predictive check to test replication for ANOVA models. }
\usage{
prior.predictive.check(n,posterior,statistic,obs=TRUE,F_n,repdata=FALSE,
                       Amat=0L,exact=0L,difmin=0L,effectsize=FALSE,seed=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n}{vector with the sample size per group (i.e., n_{jr}) for new study (i.e., y_r).}
\item{posterior}{a matrix (e.g., the output of Gibbs.ANOVA) with samples from the posterior based on the original data (i.e., y_o).}
\item{statistic}{the type of hypothesis to be calculated: "ineq" for inequality constrained means, "dif" for inequality constraints plus minimum differences between means, "exact" for specific values for the means.}
\item{obs}{logic; If FALSE, the prior predictive check does not calculate a p-value, because no observed statistic is provided. Used by the sample.size.calculator function.}
\item{F_n}{The Fbar value for the new data.}
\item{repdata}{logic; If TRUE, it uses the replication data produced by a previous run of the prior.predictive.check function.}
\item{Amat}{a p by q matrix, where p is the number of means in the ANOVA model, and q is the number of constraints to be imposed on the model. Each row represents one constraint where the parameter with the lower value according to the constraint receives the value -1, and the parameter with the higher value according to the constraint receives the value 1. Other parameters within the same row obtain the value 0.}
\item{exact}{a vector of length p, where p is the number of means in the ANOVA model, with the exact values of the constrained hypothesis.}
\item{difmin}{a vector of length q with the minimum difference per constraint as specified in \code{Amat}.}
\item{effectsize}{logical; If TRUE the values in \code{difmin} are interpreted as Cohen's d.}
\item{seed}{integer; seed value. If seed==0, no seed is set.}
}
\value{
Generates a histogram of F_sim in which F_n is indicated with a vertical line. The proportion of F_sim at the right of this line constitutes the prior predictive p-value.
\item{sumFdist}{a summary of F_sim}
\item{ppp}{the prior predictive p-value}
\item{F_sim}{a vector with F-bar values for all simulated datasets}
}
\references{
Zondervan-Zwijnenburg, M.A.J., Van de Schoot, R., & Hoijtink, H. (2017). Testing ANOVA replication by means of the prior predictive p-value. }
\author{
M. A. J. Zondervan-Zwijnenburg}
\seealso{
See also \code{\link{Gibbs.ANOVA}}, \code{\link{Fbar.ineq}}, \code{\link{Fbar.dif}}, and \code{\link{Fbar.exact}}, \code{\link{sample.size.calc}}.
}
\examples{
#analysis original data
data_o <- data.frame(y=ChickWeight$weight,g=ChickWeight$Diet)
Gibbs.ANOVA(data_o,it=75,burnin=50) #we advise >1000 iterations, >500 burnin for reliable results

#analysis new data
data_r <- data.frame(y=rnorm(660, mean(data_o$y), sd=sd(data_o$y)),g=round(runif(660,1,4)))
n.r = as.numeric(table(data_r$g))
Amat <- (rbind(c(-1,0,0,1),c(0,-1,0,1),c(0,0,-1,1)))
difmin=c(0.8,0.5,0.2)
r.F.dif.efsz <- Fbar.dif(data_r,Amat,difmin,effectsize=TRUE)

#prior predictive check
result <- prior.predictive.check(n=n.r,posterior=output_m,F_n=r.F.dif.efsz,statistic="dif",
effectsize=TRUE,Amat=Amat,difmin=difmin,repdata=FALSE,seed=1)

result$sumFdist #summary of the f(F_y_sim)
result$ppp      #the prior predictive p-value
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{nonparametric }% use one of  RShowDoc("KEYWORDS")
