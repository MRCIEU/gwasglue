% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cojo.r
\name{cojo_cond}
\alias{cojo_cond}
\title{Perform conditional analysis using GCTA COJO}
\usage{
cojo_cond(
  vcffile,
  bfile,
  snplist,
  pop,
  gcta = genetics.binaRies::get_gcta_binary(),
  workdir = tempdir(),
  threads = 1
)
}
\arguments{
\item{vcffile}{Path to vcffile}

\item{bfile}{LD reference panel}

\item{snplist}{List of rsids}

\item{pop}{EUR, ASN or AFR}

\item{gcta}{Path to gcta binary. For convenience can use default=genetics.binaRies::get_gcta_binary()}

\item{workdir}{Location to store temporary files. Default=tempdir()}

\item{threads}{Number of parallel threads. Default=1}
}
\value{
List of independent summary stats
}
\description{
For a list of fine-mapped rsids, will assign to regions and generate conditionally independent summary stats for each rsid
}
