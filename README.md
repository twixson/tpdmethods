
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tpdmethods

<!-- badges: start -->

<!-- badges: end -->

The goal of tpdmethods is to make TPD-based methods accessible to more
researchers.

Recent work in multivariate extremes is built on the so-called tail
pairwise dependence (TPD). This summary metric of dependence in the
pairwise tail shares properties with covariance and this recent work has
used these properties to construct extremal analogues to classical
statistical tools including PCA and ARMA time series models. We
implement these extremal analysis pipelines in this R package
‘tpdmethods’. The package includes datasets and vignettes for
precipitation and fire weather analyses.

## Installation

You can install the development version of tpdmethods from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak") # pak is a method to download packages from a variety of sources.
pak::pak("twixson/tpdmethods")
```

## Example

This is a basic example which shows a simple TLETS workflow. We will
perform a simplified version of the analysis in Wixson and Cooley
(2023). The scientific question is how much more likely is a fire season
like 2020 under present climate than under past climate. Wixson and
Cooley (2023) used Fire Weather Index (FWI) data to answer this
question. Their marginally-transformed data lives in the package as
`fire_weather_present`, and `fire_weather_past`. Each dataset is a
matrix with 153 rows for the days in the fire season and 20 columns for
the years.

The first step is to estimate the TPD of the time series data. Since
these data have plausibly iid replicates (each season) and are already
on Fr'echet(2) margins we will use `matrix_as_seasons = TRUE` and
`trans_marginal = FALSE` in our `tpd()` call.

``` r
library(tpdmethods)

plot(fire_weather_present[,15], 
     type = "l",
     main = "2020 Season FWI on RV Margins")
```

<img src="man/figures/README-estimate tpd-1.png" width="100%" />

``` r

fw_tpd_present <- tpd(fire_weather_present, 
                      max_lag = 30, 
                      trans_marginal = FALSE, 
                      matrix_as_seasons = TRUE)
#> [1] "Assuming this is a time series with 20 seasons."

plot(0:30, fw_tpd_present, type = "h", ylim = c(0, 1))
```

<img src="man/figures/README-estimate tpd-2.png" width="100%" />

After estimating the TPD we need to fit a model to the data. We will use
the extremes analogue to the innovations algorithm to fit a
transformed-linear moving average (TL-MA) model up to order `q = 20`. We
will visually compare the plot of the fitted TPD function to the
empirical TPD function to determine which model to retain.

``` r
set.seed(1982374)
model_coefs_present <- innovations(fw_tpd_present, max_q = 20)

plot(0:30, fw_tpd_present, type = "h", ylim = c(0, 1), 
     main = "Present TPDFs", 
     xlab = "lag", 
     ylab = "TPD value")
lines(1:30 + 0.2, maq_tpdf(model_coefs_present[[1]][15,], max_lag = 30), 
      type = "h", col = "4")
legend("topright", legend = c("empirical", "MA(15)"), 
       text.col = c(1,4), bty = "n")
abline(h = 0.1, lty = 2, col = "grey")
```

<img src="man/figures/README-fit MAs-1.png" width="100%" />

Here we can see the TL-MA(15) model fits well up to lag-15. TL-MA(q)
models have non-zero TPDFs only up to lag-q just like classical MA(q)
models. Included on the plot is a horizontal dashed line at 0.1. This is
included because of the known bias in the TPD estimator. Wixson and
Cooley (2023) considered TPD values that were consistently around this
line to be that known bias and thus continued their analysis with the
TL-MA(15).

We will repeat the previous steps with the past-climate data (from
1958-1979).

``` r
fw_tpd_past <- tpd(fire_weather_past, 
                   max_lag = 30, 
                   trans_marginal = FALSE, 
                   matrix_as_seasons = TRUE)
#> [1] "Assuming this is a time series with 20 seasons."

model_coefs_past <- innovations(fw_tpd_past, max_q = 20)

plot(0:30, fw_tpd_past, type = "h", ylim = c(0, 1), 
     main = "Past TPDFs", 
     xlab = "lag", 
     ylab = "TPD value")
lines(1:30 + 0.2, maq_tpdf(model_coefs_past[[1]][15,], max_lag = 30), 
      type = "h", col = "4")
legend("topright", legend = c("empirical", "MA(15)"), 
       text.col = c(1,4), bty = "n")
abline(h = 0.1, lty = 2, col = "grey")
```

<img src="man/figures/README-past climate-1.png" width="100%" />

Wixson and Cooley (2023) considered a high quantile of the present
climate FWI and considered a season as extreme as 2020 if it had at
least as many days above that high quantile as 2020 did. Here we use the
marginally transformed data to perform the same analysis for one
quantile (the 0.98 quantile). This analysis differs from Wixson and
Cooley because the marginal transformation removed seasonality which
included a location-parameter based climate signal. Here we are only
looking at the dependence information.

``` r
q098 <- stats::quantile(fire_weather_present, probs = 0.98)
(count_2020 <- sum(fire_weather_present[,15] > q098))
#> [1] 14
```

We note that a generated TL-MA(q) time series has an unknown marginal
distribution. This means that we will need to marginally transform the
generated seasons back to being marginally Fr'echet(2). To do this we
generate 1000 seasons all at once and use the ECDF.

``` r
set.seed(2389098)
temp_past <- gen_maq(n = 153*1000, thetas = model_coefs_past[[1]][15,1:15])
temp_present <- gen_maq(n = 153*1000, 
                        thetas = model_coefs_present[[1]][15,1:15])
```

``` r
ecdf_past    <- ecdf(temp_past)
ecdf_present <- ecdf(temp_present)

unif_past <- ecdf_past(temp_past)
unif_past <- unif_past - 0.5*min(unif_past)
unif_present <- ecdf_present(temp_present)
unif_present <- unif_present - 0.5*min(unif_present)

frechet_past    <- evd::qfrechet(p = unif_past, shape = 2)
frechet_present <- evd::qfrechet(p = unif_present, shape = 2)

seasons_past <- matrix(frechet_past, ncol = 1000)
seasons_present <- matrix(frechet_present, ncol = 1000)
```

Finally, we count the number of seasons with at least 14 days above the
high quantile.

``` r
extreme_past    <- 
  apply(seasons_past, 2, function(x){sum(x > q098) > count_2020})
extreme_present <- 
  apply(seasons_present, 2, function(x){sum(x > q098) > count_2020})

result <- data.frame(quantile = 0.98,
                     threshold = unname(q098),
                     Num_2020 = count_2020,
                     past = sum(extreme_past),
                     present = sum(extreme_present),
                     ratio = sum(extreme_present)/sum(extreme_past))
result
#>   quantile threshold Num_2020 past present     ratio
#> 1     0.98  4.599664       14  120      88 0.7333333
```

It appears that the dependence may have been slightly stronger under the
climate of fifty years ago. That hypothesis would need more careful
testing, but now you see some of the TLETS aspects of `tpdmethods`.
