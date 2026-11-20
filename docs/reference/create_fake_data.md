# Simulate synthetic LTBI diagnostic and TB progression survival data

Generates synthetic data for Bayesian Multi-Parameter Evidence Synthesis
(MPES) linking latent tuberculosis infection (LTBI) diagnostic test data
and active TB progression time-to-event data. Supports imperfect
diagnostic tests (\`sens\`, \`spec\`), multiple diagnostic tests per
person, linked or unlinked cohorts, test-positive subgroup follow-up,
covariates, and prior/data conflict scenarios (\`p_ltbi_surv !=
p_ltbi\`).

## Usage

``` r
create_fake_data(
  N_surv = 500,
  n_survey = c(1000),
  p_ltbi = 0.3,
  p_ltbi_surv = NULL,
  sens = 0.85,
  spec = 0.95,
  surv_dist = c("gompertz", "weibull", "exponential"),
  shape = 0.1,
  rate = 0.1,
  t_censor = 10,
  cohort_type = c("all", "pos_only"),
  include_covariates = FALSE,
  beta_cov = c(0.4, 0.6),
  seed = NULL
)
```

## Arguments

- N_surv:

  Integer. Sample size of the progression survival cohort (before any
  \`pos_only\` filtering). Default is \`500\`.

- n_survey:

  Integer vector. Sample sizes of aggregate LTBI prevalence surveys.
  Default is \`c(1000)\`.

- p_ltbi:

  Numeric in \`(0, 1)\`. True baseline LTBI prevalence in the diagnostic
  survey population. Default is \`0.30\`.

- p_ltbi_surv:

  Numeric in \`(0, 1)\` or \`NULL\`. True baseline LTBI prevalence in
  the progression survival cohort. Defaults to \`p_ltbi\` (no conflict).
  Set to a different value (e.g., \`0.60\`) to simulate selection bias /
  conflict between the survey and survival modules.

- sens:

  Numeric vector in \`(0, 1\]\`. Sensitivity of each diagnostic test \`k
  = 1, ..., K\`. Default is \`0.85\`.

- spec:

  Numeric vector in \`(0, 1\]\`. Specificity of each diagnostic test \`k
  = 1, ..., K\`. Default is \`0.95\`.

- surv_dist:

  Character. Parametric survival distribution for uncured progression
  times: \`"gompertz"\`, \`"weibull"\`, or \`"exponential"\`.

- shape:

  Numeric. Shape parameter for \`"gompertz"\` (can be real) or
  \`"weibull"\` (must be positive). Ignored for \`"exponential"\`.
  Default \`0.1\`.

- rate:

  Numeric (\`\> 0\`). Hazard rate / scale parameter \`lambda\`. Default
  \`0.1\`.

- t_censor:

  Numeric (\`\> 0\`). Maximum administrative follow-up / censoring time.
  Default is \`10\`.

- cohort_type:

  Character. Either \`"all"\` (all screened individuals are followed for
  active TB progression) or \`"pos_only"\` (only individuals testing
  positive on test 1 are followed for progression, so false positives
  are cured and the effective uncured fraction is the positive
  predictive value, PPV).

- include_covariates:

  Logical. If \`TRUE\`, simulates two baseline covariates (\`age_z\` and
  \`high_inc_origin\`) affecting \`logit(p_ltbi)\`.

- beta_cov:

  Numeric vector of length 2 used when \`include_covariates = TRUE\`.
  Default is \`c(0.4, 0.6)\`.

- seed:

  Optional integer random seed for reproducibility.

## Value

A list with components:

- \`progression_data\`:

  Tibble of individual-level survival follow-up (\`t\`, \`d\`,
  \`z_true\`, test results \`x1\`, ..., and optional covariates).

- \`prevalence_data\`:

  Tibble of aggregate binomial survey data (\`y\`, \`n\`, \`test_id\`)
  plus normal logit plug-in summaries (\`mu_hat\`, \`sigma_hat\`).

- \`diagnostic_ipd\`:

  Tibble of individual-level test results (\`x1\`, ...) for the
  progression cohort.

- \`truth\`:

  List of true data-generating parameters (\`p_ltbi\`, \`p_ltbi_surv\`,
  \`sens\`, \`spec\`, \`ppv_surv\`, \`surv_dist\`, \`shape\`, \`rate\`,
  \`cohort_type\`).
