
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EvidSynthTB

<!-- badges: start -->

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

## Model

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

![](docs/evid-synthesis-method/DAG-full_model.PNG)

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(EvidSynthTB)
## basic example code
```
