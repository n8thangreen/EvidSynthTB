# Changelog

## EvidSynthTB 0.2.0

- Migrated Bayesian backend to `cmdstanr` and `posterior` with
  modernized Stan models in `inst/stan/` (`ltbi_module.stan` and
  `evidsynth_modular.stan`).
- [`compare_evidsynth()`](https://n8thangreen.github.io/EvidSynthTB/reference/compare_evidsynth.md),
  [`plot_model_comparison()`](https://n8thangreen.github.io/EvidSynthTB/reference/plot_model_comparison.md),
  and
  [`plot_survival_comparison()`](https://n8thangreen.github.io/EvidSynthTB/reference/plot_survival_comparison.md)
  added to fit and compare multiple synthesis strategies side-by-side.
- [`create_fake_data()`](https://n8thangreen.github.io/EvidSynthTB/reference/create_fake_data.md)
  implemented to simulate linked or unlinked LTBI diagnostic and active
  TB progression survival cohorts with imperfect tests (`sens`, `spec`),
  multiple tests per person, `pos_only` subgroup follow-up, covariates,
  and prior/data conflict scenarios.
- [`diagnose_feedback()`](https://n8thangreen.github.io/EvidSynthTB/reference/diagnose_feedback.md)
  and
  [`evaluate_comparison()`](https://n8thangreen.github.io/EvidSynthTB/reference/evaluate_comparison.md)
  added to quantify inter-module conflict, feedback shift relative to
  the cut model, survival WAIC, parameter bias/RMSE/coverage, and
  survival curve integrated absolute error.
- [`evidsynth_fit()`](https://n8thangreen.github.io/EvidSynthTB/reference/evidsynth_fit.md)
  expanded to support six modular synthesis strategies (`"joint"`,
  `"cut"`, `"independent"`, `"melding"`, `"elastic"`, `"tempered"`),
  both `"population"` and `"latent_class"` coupling, `"all"` and
  `"pos_only"` progression cohorts, and `"gompertz"`, `"weibull"`, and
  `"exponential"` survival distributions.
- [`load_abubakar_2018()`](https://n8thangreen.github.io/EvidSynthTB/reference/load_abubakar_2018.md)
  and
  [`load_predict_migrant()`](https://n8thangreen.github.io/EvidSynthTB/reference/load_predict_migrant.md)
  added to load the reconstructed UK PREDICT-TB Kaplan-Meier dataset and
  cleaned migrant cohort dataset.
- [`plot_progression()`](https://n8thangreen.github.io/EvidSynthTB/reference/plot_progression.md)
  updated to support both uncured (`"uncured"`) and population mixture
  (`"population"`) survival curves from
  [`posterior::draws_df`](https://mc-stan.org/posterior/reference/draws_df.html)
  output.
