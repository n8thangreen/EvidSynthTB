# Fit Bayesian Multi-Parameter Evidence Synthesis Model for LTBI and TB Progression

Links latent tuberculosis infection (LTBI) diagnostic test / survey data
with active TB progression data using a parametric mixture cure model in
CmdStan. Supports both \*\*published aggregate / surveillance report
tables\*\* (cumulative subgroup counts, UKHSA/ETS annual interval
counts, and digitized Kaplan-Meier life tables) and \*\*individual-level
survival / linked test data\*\*, with imperfect diagnostic tests
(\`sens\`, \`spec\`) and six modular synthesis strategies (\`"joint"\`,
\`"cut"\`, \`"independent"\`, \`"melding"\`, \`"elastic"\`,
\`"tempered"\`).

## Usage

``` r
evidsynth_fit(
  prevalence_data = NULL,
  progression_data,
  model = c("cut", "independent", "joint", "melding", "elastic", "tempered"),
  coupling = c("population", "latent_class"),
  cohort_type = c("all", "pos_only"),
  surv_dist = c("gompertz", "weibull", "exponential"),
  imperfect_test = FALSE,
  covariates = NULL,
  t_horizon = 5,
  priors = list(),
  tau_couple = 0.5,
  w_temper = 0.5,
  t_lim = 20,
  chains = 2,
  n_iter = 1000,
  n_burnin = 500,
  n_thin = 1,
  n_cut_subsamples = 25,
  seed = 123,
  refresh = 0
)
```

## Arguments

- prevalence_data:

  Data frame or tibble of aggregate LTBI survey data (columns \`y\` or
  \`ltbi\`/\`n_pos\` and \`n\` or \`pop\`/\`n_migrant\`, plus optional
  \`test_id\` and covariates) and/or normal logit plug-in summaries
  (\`mu_hat\` and \`sigma_hat\`). Can be \`NULL\` when
  \`progression_data\` already contains group-level LTBI counts (\`y\` /
  \`ltbi\` / \`n_pos\`) or individual-level test results (\`x\` /
  \`x1\`, \`x2\`, ...).

- progression_data:

  Data frame or tibble of active TB progression data in one of four
  formats:

  1.  \*\*Published cumulative subgroup summary table:\*\* columns \`n\`
      (or \`pop\` / \`n_migrant\`), \`d\` (or \`tb\` / \`n_tb\`),
      optional \`y\` (or \`ltbi\` / \`n_pos\`), optional \`t_end\`
      (follow-up duration, defaulting to \`t_horizon\`), and optional
      subgroup covariates.

  2.  \*\*UKHSA / ETS annual or interval surveillance table:\*\* columns
      \`t_start\` (\`t0\`), \`t_end\` (\`t1\`), \`d\` (\`tb\` / \`n_tb\`
      / \`n.event\`), and \`n\` (\`pop\` / \`n_migrant\`, initial cohort
      denominator).

  3.  \*\*Published Kaplan-Meier life table:\*\* columns
      \`t_start\`/\`t_end\` (or \`time\`), \`n_risk\` (\`n.risk\`), and
      \`d\` (\`n.event\`).

  4.  \*\*Individual-level survival data (IPD):\*\* columns \`t\`
      (follow-up time \`\> 0\`) and \`d\` (\`1\` = active TB
      progression, \`0\` = right-censored), plus optional individual
      diagnostic test columns (\`x\` or \`x1\`, \`x2\`, ...).

- model:

  Character. Synthesis model strategy:

  \`"joint"\`

  :   Full joint synthesis sharing baseline LTBI prevalence
      (\`alpha_surv = alpha_diag\`).

  \`"cut"\`

  :   True two-stage Bayesian cut model (Plummer 2015): fits Stage 1 on
      the LTBI diagnostic/survey data alone, then propagates posterior
      draws of prevalence and test accuracy into Stage 2 survival
      inference with zero feedback from survival to prevalence.

  \`"independent"\`

  :   Disjoint modules with separate independent prevalence parameters
      \`prev_diag\` and \`prev_cf\`.

  \`"melding"\`

  :   Two-stage normal portal approximation (\`mu_hat\`, \`sigma_hat\`
      on \`logit(pi)\`) combined with the progression survival
      likelihood.

  \`"elastic"\`

  :   Hierarchical coupling where \`logit(prev_cf) ~
      Normal(logit(prev_diag), tau_couple)\`.

  \`"tempered"\`

  :   Power-likelihood synthesis weighting the progression survival
      likelihood by \`w_temper\` in \`(0, 1\]\`.

- coupling:

  Character. How diagnostic test results and progression times are
  coupled when individual-level linked test data are provided:
  \`"population"\` (modular population-level prevalence link) or
  \`"latent_class"\` (individual-level marginalization over each
  person's latent LTBI state \`Z_i\`).

- cohort_type:

  Character. \`"all"\` if all screened individuals are followed for
  progression (uncured fraction is \`prev_cf\`), or \`"pos_only"\` if
  only test-positive individuals are followed (uncured fraction is the
  positive predictive value, \`ppv_cf\`).

- surv_dist:

  Character. Parametric uncured survival family: \`"gompertz"\`,
  \`"weibull"\`, or \`"exponential"\`.

- imperfect_test:

  Logical. If \`TRUE\`, estimates test sensitivity (\`sens\`) and
  specificity (\`spec\`) with \`Beta(a_sens, b_sens)\` and
  \`Beta(a_spec, b_spec)\` priors. If \`FALSE\`, fixes them to
  \`sens_fixed\` and \`spec_fixed\`.

- covariates:

  Optional character vector of covariate column names to include in the
  linear predictor for \`logit(pi)\`.

- t_horizon:

  Numeric (\`\> 0\`). Default follow-up duration (in years) for
  published cumulative subgroup summary tables when \`t_end\` is not
  explicitly provided as a column. Default is \`5\`.

- priors:

  Named list of prior hyperparameters overriding defaults (\`mu_alpha\`,
  \`sigma_alpha\`, \`sigma_beta\`, \`mu_shape\`, \`sigma_shape\`,
  \`a_lambda\`, \`b_lambda\`, \`a_sens\`, \`b_sens\`, \`a_spec\`,
  \`b_spec\`, \`sens_fixed\`, \`spec_fixed\`).

- tau_couple:

  Numeric (\`\> 0\`). Coupling standard deviation on the logit scale
  when \`model = "elastic"\`. Default is \`0.5\`.

- w_temper:

  Numeric in \`\[0, 1\]\`. Power-likelihood weight on the survival
  module when \`model = "tempered"\`. Default is \`0.5\`.

- t_lim:

  Integer (\`\>= 1\`). Maximum integer time horizon for posterior
  predicted survival curves \`S_pred\[1:t_lim\]\`. Default is \`20\`.

- chains:

  Integer. Number of MCMC chains. Default is \`2\`.

- n_iter:

  Integer. Total iterations per chain (or post-warmup draws +
  \`n_burnin\`). Default is \`1000\`.

- n_burnin:

  Integer. Warmup iterations per chain. Default is \`500\`.

- n_thin:

  Integer. Thinning interval. Default is \`1\`.

- n_cut_subsamples:

  Integer. Number of Stage 1 posterior draws to propagate in Stage 2
  conditional HMC when \`model = "cut"\`. Default is \`25\`.

- seed:

  Optional integer random seed.

- refresh:

  Integer. CmdStan output refresh rate (\`0\` suppresses progress
  printing). Default is \`0\`.

## Value

An S3 list object of class \`"evidsynth_fit"\` containing:

- \`draws\`:

  A \`posterior::draws_df\` data frame of posterior samples for
  \`prev_cf\`, \`prev_diag\`, \`ppv_cf\`, \`sens\`, \`spec\`,
  \`lambda\`, \`shape\`, \`S_pred\[1:t_lim\]\`,
  \`S_pop_pred\[1:t_lim\]\`, and \`log_lik_surv\`.

- \`summary\`:

  Posterior summary table produced by \`posterior::summarise_draws()\`.

- \`fit\`:

  Underlying \`CmdStanMCMC\` object (or Stage 1/2 fit objects for
  \`"cut"\`).

- \`data\`:

  Constructed Stan data list.

- \`model\`:

  Character name of the synthesis strategy.

- \`surv_dist\`:

  Character name of the uncured survival distribution.
