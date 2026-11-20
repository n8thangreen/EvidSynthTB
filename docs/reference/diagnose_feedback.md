# Diagnose inter-module conflict, feedback shift, and survival WAIC

Quantifies information flow and potential conflict between the LTBI
diagnostic module and the active TB progression survival module across
fitted synthesis strategies, without requiring known ground-truth
parameters.

## Usage

``` r
diagnose_feedback(comparison_obj)
```

## Arguments

- comparison_obj:

  An \`"evidsynth_comparison"\` object returned by
  \[compare_evidsynth()\].

## Value

A tibble with one row per fitted model containing:

- \`model\`:

  Synthesis strategy name.

- \`prev_diag_mean\`, \`prev_diag_sd\`:

  Posterior mean and SD of LTBI prevalence in the diagnostic module.

- \`prev_cf_mean\`, \`prev_cf_sd\`:

  Posterior mean and SD of the uncured LTBI fraction in the survival
  module.

- \`feedback_shift\`:

  Shift in \`prev_diag_mean\` relative to the \`"cut"\` model (or
  \`"independent"\` if \`"cut"\` was not fitted), measuring how much
  survival data pull the LTBI prevalence estimate.

- \`conflict_diff_median\`, \`conflict_prob_pos\`:

  Posterior median of \`prev_cf - prev_diag\` and posterior probability
  \`P(prev_cf \> prev_diag \| Y)\` (informative in \`"independent"\` and
  \`"elastic"\` models).

- \`waic_surv\`, \`elpd_waic\`, \`p_waic\`:

  Widely Applicable Information Criterion (WAIC) for the progression
  survival module computed from pointwise \`log_lik_surv\`.
