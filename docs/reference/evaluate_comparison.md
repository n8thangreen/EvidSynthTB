# Evaluate parameter recovery and survival curve error against ground truth

Computes simulation recovery metrics (Bias, RMSE, Posterior SD, 95
Interval Coverage, and 95 Absolute Error (IAE) / Integrated RMSE (IRMSE)
for the uncured survival curve across models in an
\`"evidsynth_comparison"\` object.

## Usage

``` r
evaluate_comparison(comparison_obj, truth)
```

## Arguments

- comparison_obj:

  An \`"evidsynth_comparison"\` object returned by
  \[compare_evidsynth()\].

- truth:

  Either the \`truth\` list returned by \[create_fake_data()\] or a
  named numeric vector/list of true parameter values (e.g.,
  \`list(prev_diag = 0.25, prev_cf = 0.25, lambda = 0.1, shape =
  0.1)\`).

## Value

A list with two tibbles:

- \`parameter_metrics\`:

  Tibble with columns \`model\`, \`variable\`, \`truth\`, \`post_mean\`,
  \`post_median\`, \`post_sd\`, \`bias\`, \`rmse\`, \`cov_95\`,
  \`width_95\`.

- \`survival_metrics\`:

  Tibble with columns \`model\`, \`iae_surv\`, and \`irmse_surv\`
  comparing posterior median uncured survival against the true survival
  curve (if \`rate\` / \`lambda\` and \`shape\` are in \`truth\`).
