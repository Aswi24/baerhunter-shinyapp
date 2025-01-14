% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_bam_processing.R
\name{peak_analysis}
\alias{peak_analysis}
\title{Peak checking for the second coverage threshold and width.}
\usage{
peak_analysis(View_line, high_cutoff, peak_width)
}
\arguments{
\item{View_line}{A line from a RleViews object.}

\item{high_cutoff}{An integer indicating the high coverage threshold value.}

\item{peak_width}{An integer indicating the minimum sRNA length (peak width).}
}
\value{
Returns a RleViews line if it satisfies conditions.
}
\description{
This is a helper function that is used to examine if the peak had a continuous stretch of a given width that has coverage above the high cut-off value.
}
