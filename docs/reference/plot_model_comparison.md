# Plot posterior parameter comparisons across synthesis models

Creates a faceted density plot comparing the posterior distributions of
key parameters (\`prev_cf\`, \`prev_diag\`, \`lambda\`, \`shape\`,
\`sens\[1\]\`, \`spec\[1\]\`) across synthesis strategies.

## Usage

``` r
plot_model_comparison(
  comparison_obj,
  variables = c("prev_cf", "prev_diag", "lambda", "shape"),
  truth = NULL
)
```

## Arguments

- comparison_obj:

  An \`"evidsynth_comparison"\` object returned by
  \[compare_evidsynth()\].

- variables:

  Character vector of parameter names to display. Default is
  \`c("prev_cf", "prev_diag", "lambda", "shape")\`.

- truth:

  Optional named numeric vector or list of true parameter values (e.g.
  \`c(prev_cf = 0.3, prev_diag = 0.3, lambda = 0.1, shape = 0.1)\`) to
  overlay as vertical dashed reference lines.

## Value

A \`ggplot2\` object.
