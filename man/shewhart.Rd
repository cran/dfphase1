\name{shewhart}
\alias{shewhart}
\alias{shewhart.normal.limits}
\title{
  Univariate Shewhart-type control charts
}
\description{
  \code{shewhart} computes, and, optionally,  plots, 
  Shewhart-type Phase I control charts for detecting
  changes in location and scale of univariate subgrouped data.
  

  \code{shewhart.normal.limits} pre-computes
  the corresponding control limits when the in-control distribution is normal.
}
\usage{
shewhart(x, subset, 
         stat = c("XbarS", "Xbar", "S", 
                  "Rank", "lRank", "sRank",
                  "Lepage", "Cucconi"),
         aggregation = c("mean", "median"), 
         plot = TRUE, 
         FAP = 0.05,
         seed = 11642257, 
         L = 1000, 
         limits = NA)

shewhart.normal.limits(n, m, 
                       stat = c("XbarS", "Xbar", "S", 
                                "Rank", "lRank", "sRank", 
                                "Lepage", "Cucconi"),
                       aggregation = c("mean", "median"), 
                       FAP = 0.05,
                       seed = 11642257, 
                       L = 100000)
}
\arguments{
  \item{x}{
    a nxm data numeric matrix
    (n observations gathered at m time points).
  }
  \item{subset}{
    an optional vector specifying a subset of subgroups/time points to be used
  }
  \item{stat}{
    character: the control statistic[s] to use; see Details.
}
\item{aggregation}{
  character: 
  it specify how to aggregate the subgroup means and standard deviations.
  Used only when \code{stat} is \code{XbarS}, \code{Xbar} or \code{S}.
}
  \item{plot}{
        logical; if \code{TRUE}, control statistic[s] is[are] displayed. 
}
  \item{FAP}{
    numeric (between 0 and 1): desired false alarm probability.
    Unused by \code{shewhart} when \code{limits} is not \code{NA}.
}
  \item{seed}{
    positive integer; if not \code{NA}, the RNG's state is resetted
    using \code{seed}. The current \code{.Random.seed} will be
    preserved.
    Unused by \code{shewhart} when \code{limits} is not \code{NA}.
}
  \item{L}{
    positive integer: number of random permutations used to
    compute the control limits.  Unused by \code{shewhart} when \code{limits} is not \code{NA}.
}
  \item{limits}{
    numeric: a precomputed vector of control limits.
    The vector should contain \eqn{(A,B_1,B_2)}{(A,B1,B2)}
    when \code{stat=XbarS}, \eqn{(A)} when \code{stat=Xbar},
    \eqn{(B_1,B_2)}{(B1,B2)} when \code{stat=S},
    \eqn{(C,D)} when \code{stat=Rank}, \eqn{(C)} when
    \code{stat=lRank}, \eqn{(D)} when \code{stat=sRank},
    and \eqn{(E)} when \code{stat=Lepage} or \code{stat=Cucconi}.
    See Details for the definition of the critical values
    \eqn{A}, \eqn{B_1}{B1},  \eqn{B_2}{B2}, \eqn{C}, \eqn{D}
   and \eqn{E}. 
}
\item{n}{
    integer: size of each subgroup (number of observations gathered at
    each time point).
  }
\item{m}{
    integer: number of subgroups (time points).
  }
}
\details{
  The implemented control charts are:
  \itemize{
    \item{\code{XbarS}}: combination of the \code{Xbar}
    and \code{S} control charts described in the following.
    \item{\code{Xbar}}: chart based on plotting the subgroup means with control limits
    \deqn{\hat{\mu}\pm A\frac{\hat{\sigma}}{\sqrt{n}}}{mu.hat +/- A
      sigma.hat/sqrt(n)}
    where \eqn{\hat{\mu}}{mu-hat} (\eqn{\hat{\sigma}}{sigma.hat})
    denotes the estimate of the in-control mean (standard deviation)
    computed as the mean or median of the subgroup means (standard
    deviations).
    \item{\code{S}}:  chart based on plotting the (unbiased) subgroup standard deviations
    with lower control limit \eqn{B_1\hat{\sigma}}{B1 x sigma.hat}  and
    upper control limit \eqn{B_2\hat{\sigma}}{B2 x sigma.hat}.
    \item{\code{Rank}}: combination of the \code{lRank}
    and \code{sRank} control charts described in the following.
    \item{\code{lRank}}: control chart based on the standardized
    rank-sum control statistic suggested by
    Jones-Farmer \emph{et al.} (2009) for detecting changes in the location parameter.
    Control limits are of the type \eqn{\pm C}{+/- C}.
    \item{\code{sRank}}: chart based on the standardized
    rank-sum control statistic suggested by
    Jones-Farmer and Champ (2010) for detecting changes in the scale parameter.
    Control limits are of the type \eqn{\pm D}{+/- D}.
    \item{\code{Lepage}}: chart based on the Lepage control statistic
    suggested by Li \emph{et al.} (2019) for detecting changes in
    location and/or scale. There is only a upper control limit equal to \eqn{E}.
    \item{\code{Cucconi}}: chart based on the Cucconi control statistic
    suggested by Li \emph{et al.} (2020) for detecting changes in
    location and/or scale. There is only a upper control limit equal to \eqn{E}.
  }
}


\value{
  \code{shewhart} returns an invisible list with elements
  \item{\code{Xbar} }{subgroup means; this element is present only if
    \code{stat} is \code{XbarS} or \code{Xbar}.}
  \item{\code{S} }{subgroup standard deviation; this element is present only if
    \code{stat} is \code{XbarS} or \code{S}.}
  \item{\code{lRank} }{rank-based control statistics for detecting
  changes in location; this element is present only if
    \code{stat} is \code{Rank} or \code{lRank}.}
  \item{\code{sRank} }{rank-based control-statistics for detecting
  changes in scale; this element is present only if
    \code{stat} is \code{Rank} or \code{sRank}.}
  \item{\code{Lepage}, \code{W2}, \code{AB2} }{Lepage, squared Wilcoxon
      and squared Ansari-Bradley statistics; these elements are present 
      only if \code{stat} is \code{Lepage}.}
  \item{\code{Cucconi}, \code{lCucconi}, \code{sCucconi}}{Cucconi control statistic and 
      its location and scale components; 
    these elements are present only if \code{stat} is \code{Cucconi}.}   
  \item{\code{limits} }{control limits.}
  \item{\code{center}, \code{scale} }{estimates
  \eqn{\hat{\mu}}{mu.hat} and \eqn{\hat{\sigma}}{sigma.hat} of the in-control
  mean and standard deviation; these elements are present only if
  \code{stat} is \code{XbarS}, \code{Xbar} and \code{S}.}
  \item{\code{stat}, \code{L}, \code{aggregation}, \code{FAP},
    \code{seed} }{ input arguments.}

  \code{shewhart.normal.limits} returns a numeric vector
  containing the limits.
}
\note{
  \enumerate{
    \item
    If argument \code{limits} is \code{NA}, \code{shewhart}
    computes the control limits by permutation. 
    The resulting control chart are distribution-free.
    \item
    Pre-computed limits, such as those computed using
    \code{shewhart.normal.limits},  are not recommended 
    when \code{stat} is \code{XbarS}, \code{Xbar} or \code{S}.
    Indeed, the resulting control chart will not be distribution-free.
    \item
    When \code{stat} is \code{Rank}, \code{lRank},
    \code{sRank}, \code{Lepage} or \code{Cucconi} the control limits computed by
    \code{shewhart.normal.limits} are distribution-free in the class
    of all univariate continuous distributions.
    So, if user plan to apply rank-based control charts on a repeated
    number of samples of the same size,  pre-computing the control limits using
    \code{mshewhart.normal.limits}
    can reduce the overall computing time.
    }
}

\author{
\packageAuthor{dfphase1}.
}

\references{
    L. A. Jones-Farmer, V. Jordan, C. W. Champs (2009)
    \dQuote{Distribution-free Phase I control charts for subgroup location},
    \emph{Journal of Quality Technology}, \bold{41}, pp. 304--316,
   \doi{10.1080/00224065.2009.11917784}. 


    L. A. Jones-Farmer, C. W. Champ (2010) 
    \dQuote{A distribution-free Phase I control chart for subgroup scale}.
    \emph{Journal of Quality Technology}, \bold{42}, pp. 373--387,
    \doi{10.1080/00224065.2010.11917834}

    C. Li, A. Mukherjee, Q. Su (2019) 
    \dQuote{A distribution-free Phase I monitoring scheme for subgroup location and 
        scale based on the multi-sample Lepage statistic},
    \emph{Computers & Industrial Engineering}, \bold{129}, pp. 259--273,
     \doi{10.1016/j.cie.2019.01.013} 

    C. Li, A. Mukherjee, M. Marozzi (2020)
    \dQuote{A new distribution-free Phase-I procedure for bi-aspect
        monitoring based on the multi-sample Cucconi statistic},
    \emph{Computers & Industrial Engineering}, \bold{149}, 
     \doi{10.1016/j.cie.2020.106760}   

  D. C. Montgomery (2009) \emph{Introduction to Statistical Quality Control}, 6th edn. Wiley.

  P. Qiu (2013) \emph{Introduction to Statistical Process Control}. Chapman & Hall/CRC Press.

}

\examples{
# A simulated example
set.seed(12345)
y <- matrix(rt(100,3),5)
y[,20] <- y[,20]+3
shewhart(y)
shewhart(y, stat="Rank")
shewhart(y, stat="Lepage")
shewhart(y, stat="Cucconi")
# Reproduction of the control chart shown
# by Jones-Farmer et. al. (2009)
data(colonscopy)
u <- shewhart.normal.limits(NROW(colonscopy),NCOL(colonscopy), 
                            stat="lRank", FAP=0.1, L=10000)
# In Jones-Farmer et al. (2009) is estimated as 2.748
u
shewhart(colonscopy,stat="lRank",limits=u)
# Examples of control limits for comparisons
# with Li et al. (2019) and (2020) but
# using a limited number of Monte Carlo
# replications
# Lepage: in Li et al. (2019) is estimated as 11.539
shewhart.normal.limits(5, 25, stat="Lepage", L=10000)
# Cucconi: in Li et al. (2020) is estimated as 0.266
shewhart.normal.limits(5, 25, stat="Cucconi", L=10000)
}

\keyword{htest}
\keyword{statistical process monitoring}
\keyword{statistical process control}
\keyword{distribution-free methods}
