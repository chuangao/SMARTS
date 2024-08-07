% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SMARTS.R
\name{get_count_per_bin}
\alias{get_count_per_bin}
\title{cut switching time into bins, and record the numbers}
\usage{
get_count_per_bin(swi, cont, breaks, swi_time = swi_time)
}
\arguments{
\item{breaks}{a vector of quantiles used to break switching times into quantiles}
}
\value{
a list containing the switchers, continuers, the cut of the bins, number of switchers and continuers in each bin
}
\description{
cut switching time into bins, and record the numbers
}
