% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SMARTS.R
\name{random_assign}
\alias{random_assign}
\title{Randomly draw switching time from switchers and assign to continuers}
\usage{
random_assign(
  data,
  nbin = 10,
  seed = 123,
  swi_time = "xoyrs",
  cens_time = "censyrs"
)
}
\arguments{
\item{data}{a list, data$cont is the continuers, data$swi is the switchers}

\item{nbin}{number of bins to divide the switchers into}

\item{swi_time}{the switching time, pre-calculated as the time length from treatment initiation to switching. For the input data, this column only exist in the switchers, a new column with the same name will be created for the continuers after assignment}
}
\value{
a list of assigned and unassigned continuer and switchers
\itemize{
\item assigned: the continuers that are assigned switching time successfully, and the switchers that are used to assigned the continuers
\item unassigned: the continuers that are assigned switching times unsuccessfully, and the switchers that are never used 
}
}
\description{
Randomly draw switching time from switchers and assign to continuers
}
\examples{
data_tmp <- sim_switch(seed = 123, n = 5000)
data <- list(cont = data_tmp \%>\% filter(swi == 0),
             swi = data_tmp \%>\% filter(swi == 1))
data_rand <- random_assign(data, nbin=10, seed=123, swi_time='xoyrs', cens_time = 'censyrs')              
}
