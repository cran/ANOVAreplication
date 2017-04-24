\name{sample.size.calc}
\alias{sample.size.calc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sample size calculator for the prior.predictive.check function}
\description{
Calculates the required sample size for a new study to conduct the prior predictive check with sufficient statistical power.}
\usage{
sample.size.calc(start_n,itmax=10,nmax=600,powtarget=.825,powmargin=.025,posterior,g.m,
                 statistic,Amat=0L,exact=0L,difmin=0L,effectsize=FALSE,alpha=.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start_n}{integer; the starting value for the sample size per group.}
  \item{itmax}{integer; the maximum number of iterations for the function.}
  \item{nmax}{integer; the maximum total sample size to evaluate.}
  \item{powtarget}{integer; the target power for which the sample size is to be obtained.}
  \item{powmargin}{integer; the margin around the target power for which results are to be returned.}
  \item{posterior}{a matrix (e.g., the output of Gibbs.ANOVA) with samples from the posterior based on the original data (i.e., y_o).}
  \item{g.m}{integer; the population value for the equal means in the alternative distribution. We advice to specify the grand mean of the study variables in the original dataset.}
\item{statistic}{the type of hypothesis to be calculated: "ineq" for inequality constrained means, "dif" for inequality constraints plus minimum differences between means, "exact" for specific values for the means.}
\item{Amat}{a p by q matrix, where p is the number of means in the ANOVA model, and q is the number of constraints to be imposed on the model. Each row represents one constraint where the parameter with the lower value according to the constraint receives the value -1, and the parameter with the higher value according to the constraint receives the value 1. Other parameters within the same row obtain the value 0.}
\item{exact}{a vector of length p, where p is the number of means in the ANOVA model, with the exact values of the constrained hypothesis.}
\item{difmin}{a vector of length q with the minimum difference per constraint as specified in \code{Amat}.}
\item{effectsize}{logical; If TRUE the values in \code{difmin} are interpreted as Cohen's d.}
  \item{alpha}{integer; the level of alpha that should be taken into account while calculating the required sample size.}
}

\value{
Prints iterations while calculating. Prints a matrix with two columns. The first column contains the sample size per group and the second column the associated power. Furthermore, sample.size.calc produces a histogram with the null (i.e., the red distribution) and alternative distribution (i.e., the blue distribution) for the last iteration afterwards. The vertical line indicates rej.value (i.e., the 1-alpha'th percentile of the null distribution). The proportion of the alternative distribution on the right side of rej.value constitutes the statistical power.
%\item{rej.value}{integer; the 1-alpha'th percentile of the null distribution. The proportion of the alternative distribution larger than the rejection value contitutes the statistical power.}
%\item{Npower.l}{a matrix with two columns. The first column contains the sample size per group and the second column the associated power.}
}
\references{
Zondervan-Zwijnenburg, M.A.J., Van de Schoot, R., & Hoijtink, H. (2017). Testing ANOVA replication by means of the prior predictive p-value. }
\author{
M.A.J. Zondervan-Zwijnenburg}
\seealso{
See also \code{\link{prior.predictive.check}}, \code{\link{Fbar.ineq}}, \code{\link{Fbar.dif}}, and \code{\link{Fbar.exact}}.
}
\examples{
#analysis original data
data_o <- data.frame(y=ChickWeight$weight,g=ChickWeight$Diet)
Gibbs.ANOVA(data_o,it=75,burnin=50) #we advise >1000 iterations, >500 burnin for reliable results
Amat <- (rbind(c(-1,0,0,1),c(0,-1,0,1),c(0,0,1,-1)))
difmin=c(0.8,0.5,0.2)

#sample size calculator
sample.size.calc(start_n=60, powtarget=.825,powmargin=.025,posterior=output_m,
                 g.m=mean(data_o$y),statistic="dif",effectsize=TRUE,Amat=Amat,difmin=difmin)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{nonparametric }