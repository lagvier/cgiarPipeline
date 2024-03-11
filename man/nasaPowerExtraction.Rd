\name{nasaPowerExtraction}
\alias{nasaPowerExtraction}

\title{
Extraction of weather data
}
\description{
Extracts weather data from NASAPOWER database for a point of coordinates and a time interval.
}
\usage{
nasaPowerExtraction(LAT,LONG,
                    date_planted,date_harvest,
                    environments, temporal="hourly")
}

\arguments{
  \item{LAT}{
Latitude coordinate.
}
\item{LONG}{
Longitude coordinate.
}
\item{date_planted}{
Planting date.
}
\item{date_harvest}{
Harvesting date.
}
\item{environments}{
Environments.
}
\item{temporal}{
Time interval to use for extraction. 
}

}
\details{

This is a very simple function that extracts weather data from NASAPOWER.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{the same input object with an additional traits.}
}
}

\references{

None

}

\examples{
####=========================================####

# library(dplyr)
# final <- nasaPowerExtraction(LAT=0,LONG=1,
#                              date_planted=Sys.Date()-20,date_harvest=Sys.Date()-19,
#                              environment="A", temporal="hourly")
# head(final$WTH)
# head(final$descriptive)

}

