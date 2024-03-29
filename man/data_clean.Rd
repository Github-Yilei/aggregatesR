% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ternary_plot.R
\name{data_clean}
\alias{data_clean}
\title{data_clean}
\usage{
data_clean(
  otu,
  design,
  type = c("relative", "absolute"),
  threshold = 0.001,
  times = 100
)
}
\arguments{
\item{otu}{a tabbed text files in which OTUs are rows and samples are columns,
one entry in the table is usually a number of reads, or a frequency in the range 0.0 to 1.0.}

\item{design}{a data table in which fixing the header of sample as 'sampleID' and the header of group as 'Group', respectively.}

\item{type}{the type of otu table:  relative/absolute abundance otu-table}

\item{times}{the zoom scale of the points.}
}
\value{
cleaned data
}
\description{
buildding a stanard table for ggtern
}
\details{
clean data
}
\examples{
file_path <- system.file("extdata", package = "aggregatesR" )
otutab <- read.delim(paste(file_path, "otutab.txt", sep = "/"), header=T, row.names=1)
design <- read.delim(paste(file_path, "metadata.txt", sep = "/"), header=T, row.names=NULL)
design = design[ , c("SampleID","Group")]
otu_tern <- data_clean(otutab, design, type="absolute", threshold=0.001, times=100)
head(otu_tern,n=3)
p <- ggtern(data = otu_tern, aes(x = KO, y = OE, z = WT)) +
geom_point(aes(size = size), alpha = 0.8, show.legend = T) +
scale_size(range = c(0, 6)) + geom_mask() +
guides(colour = "none") + theme_bw() +
theme(axis.text = element_blank(), axis.ticks = element_blank())
}
\author{
Yilei
}
