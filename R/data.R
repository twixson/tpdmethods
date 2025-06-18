#' Fire Weather Index Data - Present period
#'
#' Fire Weather Index (FWI) data after marginal transformation to Fr\'echet(2)
#'    margins. These data are from the gridbox that contains the Harbison
#'    Meadow RAWS weather station near Grand Lake, Colorado. They are data from
#'    June through October of 2002-2021. Weather data were downloaded from the
#'    ERA5 data on the Copernicus Climate Change Service. A `matrix` with 153
#'    rows (days of the fire season) and 20 columns (years).
#'
#' @name fire_weather_present
#' @docType data
#' @references
#'  \insertRef{wixson2023}{tpdmethods}
#'
#' @source <https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation>
NULL

#' Fire Weather Index Data - Past period
#'
#' Fire Weather Index (FWI) data after marginal transformation to Fr\'echet(2)
#'    margins. These data are from the gridbox that contains the Harbison
#'    Meadow RAWS weather station near Grand Lake, Colorado. They are data from
#'    June through October of 1959 - 1978. Weather data were downloaded from the
#'    ERA5 data on the Copernicus Climate Change Service. A `matrix` with 153
#'    rows (days of the fire season) and 20 columns (years).
#'
#' @name fire_weather_past
#' @docType data
#' @references
#'  \insertRef{wixson2023}{tpdmethods}
#'
#' @source <https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation>
NULL


#' Financial Data from the Kenneth R. French Data Library.
#'
#' These data are the losses in 30 different sectors of the economy (e.g.,
#'    food, oil, and health) from 1970 through 2023. A `matrix` with 13599 rows
#'    (days financial data were recorded) and 31 columns (first column is date,
#'    the rest are the sectors).
#'
#' @name financial_data
#' @docType data
#' @references
#'  \insertRef{kiriliouk2022}{tpdmethods}
#'
#' @source <https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html>
NULL
