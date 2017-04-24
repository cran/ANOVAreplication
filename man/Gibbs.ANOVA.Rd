\name{Gibbs.ANOVA}
\alias{Gibbs.ANOVA}
\title{Gibbs sampler}
\description{
Samples from the posterior distribution of the data by means of a Gibbs sampler (derived from Lynch, 2007, p. 170-172).}
\usage{
Gibbs.ANOVA(data,it=5000,burnin=500)
}
\arguments{
\item{data}{a data frame with a variable y and a variable g, where y is the dependent variable, and g is the grouping variable for the ANOVA.}
\item{it}{the number of (post-burnin) iterations for each of the two chains. The default uses 5.000 iterations.}
\item{burnin}{the number of iterations for the function to use for the burnin phase in each of the two chains. The default uses 500 burnin iterations.}
}
\value{
\item{output_m}{a matrix with all samples from the posterior for each parameter.}

Returns a matrix with the mean, median, and standard deviation (in columns) for the it samples from the conditional posterior distributions of the group means and pooled standard deviation (rows).

Produces traceplots of each parameter and the associated samples from the posterior distribution.
}
\references{
Lynch, S. (2007). Introduction to applied Bayesian statistics and estimation for social scientists. New York, NY: Springer. doi: 10.1007/978-0-387-71265-9

Zondervan-Zwijnenburg, M.A.J., Van de Schoot, R., & Hoijtink, H. (2017). Testing ANOVA replication by means of the prior predictive p-value.}
\author{
M. A. J. Zondervan-Zwijnenburg
}

\examples{
data <- data.frame(y=ChickWeight$weight,g=ChickWeight$Diet)
Gibbs.ANOVA(data)
}
\keyword{univar}
