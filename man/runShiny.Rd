\name{runShiny}
\alias{runShiny}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function to launch Shiny application to test replication of ANOVA results
}
\description{
Launches a Shiny application for the replication test.
}
\usage{
runShiny()
}
\value{
In the Shiny application Gibbs.ANOVA can be ran, which prints covergence plots, samples of the posterior, and a summary table.

Subsequently, the sample size calculator can be used, which prints a matrix with two columns. The first column contains the sample size per group and the second column the associated power. Furthermore, sample.size.calc produces a histogram to illustrate power as evaluated in the last iteration. Detailed descriptions are provided in the Shiny application.

If information for a new study is provided, the prior predictive check can be used to compute the prior predictive p-value. The prior predictive check generates a histogram of F_sim in which F_n is indicated with a vertical line. The proportion of F_sim at the right of this line constitutes the prior predictive p-value. The user can download this histogram and sumFdist: a summary of F_sim.
}
\seealso{
See also \code{\link{Gibbs.ANOVA}}, \code{\link{sample.size.calc}}, \code{\link{prior.predictive.check}}, \code{\link{Fbar.ineq}}, \code{\link{Fbar.dif}}, and \code{\link{Fbar.exact}}.
}
\keyword{htest}
\keyword{nonparametric }
