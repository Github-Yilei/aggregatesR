% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ternary_plot.R
\name{enrich_data}
\alias{enrich_data}
\title{enrich_data}
\usage{
enrich_data(otu, design, p.value = 0.05, adjust.method = "fdr")
}
\arguments{
\item{otu}{a tabbed text files in which OTUs are rows and samples are columns,
one entry in the table is a number of reads.}

\item{design}{a data table in which fixing the header of sample as 'sampleID' and the header of group as 'Group', respectively.}
}
\value{
cleaned data
}
\description{
the differential analysis of otu table
}
\details{
enrich data
}
\examples{
enrich_index <- enrich_data(otutab, design, p.value = 0.05)
plot_data <- merge(otu_tern, enrich_index, by = "OTUs", all.x = T)
p <- ggtern(data = plot_data, aes(x = KO, y = OE, z = WT)) +
geom_mask() +
geom_point(aes(size=size, color=enrich),alpha=0.8) +
guides(size="none") +theme_bw() +
theme(axis.text=element_blank(),
axis.ticks=element_blank())
}
\author{
Yilei
}
