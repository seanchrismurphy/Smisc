\name{Smisc-package}
\alias{Smisc-package}
\alias{Smisc}
\docType{package}
\title{
Miscellaneous useful functions
}
\description{
Right now, this is just a collection of useful functions for me. It may be more at a later date.
}
\details{
\tabular{ll}{
Package: \tab Smisc\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-05-26\cr
License: \tab GPL-2 \cr
}
Important functions: 

Datamerge() which merges two datasets with some different columns by creating empty columns so that they'll match (note this may be redundant with plyr::join)

corstars() which creates pretty correlation tables

spss.load() is a wrapper around hadley's read_sav that saves a few lines

getdesc() gets the descriptions from each column of an spss file that has been loaded with spss.load or read_sav

qual.load() loads up .csv files downloaded from qualtrics and deals with the crappy qualtrics formatting, putting original question names in a description attribute that can be accessed with description()

}
\author{
Sean Murphy

Maintainer: Sean Murphy <seanchrismurphy@gmail.com>
}
