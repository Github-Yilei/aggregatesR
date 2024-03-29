% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_aspera_fq.R
\name{set_aspera_fq}
\alias{set_aspera_fq}
\title{set_aspera_fq}
\usage{
set_aspera_fq(path, ena_tsv, ascp, step)
}
\arguments{
\item{path}{the absolute path of the output file.}

\item{ena_tsv}{the ENA file.}

\item{ascp}{the absolute path of the ascp package.}

\item{step}{the jobs of each file.}

\item{openssh}{the openssh of ascp package(at etc file).}

\item{fastq_cmd}{the command line for downloading fastq data(please checking it by hand).}
}
\description{
setting the command lines of downloading fastq data from
ENA(\url{https://www.ebi.ac.uk/ena/browser/home}) by using
IBM Aspera.
}
\examples{
path <- 'C:/Users/Admin/Desktop/'
ena_tsv <- system.file("extdata", "filereport_read_run_PRJNA648176_tsv.txt", package = "aggregatesR" )
ascp <- '~/miniconda3/pkgs/aspera-cli-3.9.1-0'
set_aspera_fq(path, ena_tsv, ascp, 10)
}
