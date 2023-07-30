
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EvidSynthTB: Evidence synthesis of TB data

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of EvidSynthTB is to use Bayesian Multi-Parameter Evidence
Synthesis (MPES) to estimate TB epidemiological parameters LTBI
prevalence and active TB progression rate.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("n8thangreen/EvidSynthTB")
```

## Data

Input data sets used in this analysis are:

- Enhanced TB Surveillance (ETS)
- PREDICT-TB

## General model

We want to obtain posterior distributions for LTBI prevalence, `pl`, and
active TB activation rate, `lambda`. The other model parameters are:

- `Xm1`: Cohort size (observed)
- `Xp1`: Number positive test results (observed)
- `Xl1`: Number latent TB (unobserved)
- `XTB1`: Number active TB (observed)
- `Xm2`: Cohort size (observed)
- `XTB2`: Number active TB (observed)
- `p_pos`: Test positivity (functional)
- `pTB`: Probability active TB (functional)
- `sens`, `spec`: Test sensitivity and specificity (prior)

A Directed Acyclic Graph of the model is given below.

![](man/figures/DAG-full_model.PNG)

## Example

This is a basic example which shows you how to solve a common problem.
Fit markov melding model using Stan with simulated artificial LTBI and
TB progression data without covariate (age or ethnicity).

``` r
library(rstan)
library(shinystan)
library(purrr)
library(readr)
library(dplyr)
library(EvidSynthTB)
```

We generate some fake data.

``` r
# sample size
N <- 500

# disease-free censoring additional time
t_offset <- 10

# prevalence of ltbi
p_ltbi <- 0.3

shape0 <- 0.1
rate0 <- 0.1

# progression from LTBI to active tb times
progression_dat <-
  data.frame(
    t = round(flexsurv::rgompertz(N, shape = shape0, rate = rate0), 3)) |>
    mutate(x = rbinom(n = N, size = 1, prob = p_ltbi),   # ltbi status
           ## observe all progression times
           d = ifelse(x == 1, 1, 0),                     # censoring status
           # d = sample(c(0,1), size = N, replace = TRUE),
           # d = ifelse(x == 0, 0, d),
           t = ifelse(x == 0, t + t_offset, t)) |>       # time
  as_tibble()

# prevalence of LTBI
prevalence_dat <-
  data.frame(mu_hat = boot::logit(rnorm(n = 10, mean = 0.3, sd = 0.05)),
             sigma_hat = 0.01)
```

We can now use the prevalence data and progression data to obtain
posterior distributions.

``` r
out <- evidsynth_fit(prevalence_dat, progression_dat)
#> 
#> SAMPLING FOR MODEL 'stan_output_fake_markov_melding' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000482 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 4.82 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: WARNING: There aren't enough warmup iterations to fit the
#> Chain 1:          three stages of adaptation as currently configured.
#> Chain 1:          Reducing each adaptation stage to 15%/75%/10% of
#> Chain 1:          the given number of warmup iterations:
#> Chain 1:            init_buffer = 15
#> Chain 1:            adapt_window = 75
#> Chain 1:            term_buffer = 10
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  101 / 2000 [  5%]  (Sampling)
#> Chain 1: Iteration:  300 / 2000 [ 15%]  (Sampling)
#> Chain 1: Iteration:  500 / 2000 [ 25%]  (Sampling)
#> Chain 1: Iteration:  700 / 2000 [ 35%]  (Sampling)
#> Chain 1: Iteration:  900 / 2000 [ 45%]  (Sampling)
#> Chain 1: Iteration: 1100 / 2000 [ 55%]  (Sampling)
#> Chain 1: Iteration: 1300 / 2000 [ 65%]  (Sampling)
#> Chain 1: Iteration: 1500 / 2000 [ 75%]  (Sampling)
#> Chain 1: Iteration: 1700 / 2000 [ 85%]  (Sampling)
#> Chain 1: Iteration: 1900 / 2000 [ 95%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.794 seconds (Warm-up)
#> Chain 1:                34.728 seconds (Sampling)
#> Chain 1:                35.522 seconds (Total)
#> Chain 1:
```

We can view the output.

``` r
library(reshape2)
library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 4.2.3

plot_progression(out)
```

<img src="man/figures/README-plots-1.png" width="100%" />

``` r

# LTBI prevalence
stan_output <- extract(out$fit)

hist(stan_output$prev_cf, breaks = 40)
```

<img src="man/figures/README-plots-2.png" width="100%" />
