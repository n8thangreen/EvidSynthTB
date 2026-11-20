# Compare multiple Bayesian evidence synthesis models

Fits multiple modular/joint synthesis strategies (e.g., \`"joint"\`,
\`"cut"\`, \`"independent"\`, \`"melding"\`, \`"elastic"\`,
\`"tempered"\`) to the same dataset and collects tidy posterior
summaries and survival curves for side-by-side comparison.

## Usage

``` r
compare_evidsynth(
  prevalence_data = NULL,
  progression_data,
  models = c("joint", "cut", "independent", "melding", "elastic", "tempered"),
  ...
)
```

## Arguments

- prevalence_data:

  Aggregate survey or plug-in prevalence data passed to
  \[evidsynth_fit()\].

- progression_data:

  Individual-level survival (and optional linked test) data passed to
  \[evidsynth_fit()\].

- models:

  Character vector of synthesis strategies to fit. Default is
  \`c("joint", "cut", "independent", "melding", "elastic",
  "tempered")\`.

- ...:

  Additional arguments passed to \[evidsynth_fit()\] (such as
  \`imperfect_test\`, \`coupling\`, \`cohort_type\`, \`surv_dist\`,
  \`n_iter\`, \`n_burnin\`, \`chains\`, \`seed\`).

## Value

A list of class \`"evidsynth_comparison"\` with components:

- \`fits\`:

  Named list of \`"evidsynth_fit"\` objects.

- \`summary_table\`:

  Combined tibble of posterior summaries across all fitted models.

- \`draws_long\`:

  Combined long-format tibble of key parameter posterior draws across
  models.
