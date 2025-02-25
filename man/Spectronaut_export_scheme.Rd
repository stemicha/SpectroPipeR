% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpectroPipeR_SN_report_scheme.R
\name{Spectronaut_export_scheme}
\alias{Spectronaut_export_scheme}
\title{SpectroPipeR - Spectronaut report export scheme}
\usage{
Spectronaut_export_scheme(output_location = "")
}
\arguments{
\item{output_location}{folder location where the SpectroPipeR_report.rs report export scheme should be copied}
}
\value{
Generates a Spectronaut report export template that can be imported into Spectronaut, encompassing all necessary columns SpectroPipeR needs for the Spectronaut export of the analyzed data.
}
\description{
SpectroPipeR - Spectronaut report export scheme
}
\examples{

Spectronaut_export_scheme(output_location = "../SpectroPipeR_test_folder")

}
