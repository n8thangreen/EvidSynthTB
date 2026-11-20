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
  and
  [`create_fake_published_data()`](https://n8thangreen.github.io/EvidSynthTB/reference/create_fake_published_data.md)
  implemented to simulate individual-level cohorts, published subgroup
  summary tables, and UKHSA-style annual surveillance progression tables
  with imperfect tests (`sens`, `spec`), covariates, and prior/data
  conflict scenarios.
- [`diagnose_feedback()`](https://n8thangreen.github.io/EvidSynthTB/reference/diagnose_feedback.md)
  and
  [`evaluate_comparison()`](https://n8thangreen.github.io/EvidSynthTB/reference/evaluate_comparison.md)
  added to quantify inter-module conflict, feedback shift relative to
  the cut model, survival WAIC, parameter bias/RMSE/coverage, and
  survival curve integrated absolute error.
- [`evidsynth_fit()`](https://n8thangreen.github.io/EvidSynthTB/reference/evidsynth_fit.md)
  expanded to natively accept both published aggregate/interval/KM
  tables and individual-level survival data across six modular synthesis
  strategies (`"joint"`, `"cut"`, `"independent"`, `"melding"`,
  `"elastic"`, `"tempered"`), `"population"` and `"latent_class"`
  coupling, `"all"` and `"pos_only"` progression cohorts, and
  `"gompertz"`, `"weibull"`, and `"exponential"` survival distributions.
- [`load_abubakar_km()`](https://n8thangreen.github.io/EvidSynthTB/reference/load_abubakar_km.md),
  [`load_predict_aggregated()`](https://n8thangreen.github.io/EvidSynthTB/reference/load_predict_aggregated.md),
  [`load_abubakar_2018()`](https://n8thangreen.github.io/EvidSynthTB/reference/load_abubakar_2018.md),
  and
  [`load_predict_migrant()`](https://n8thangreen.github.io/EvidSynthTB/reference/load_predict_migrant.md)
  added to load published Kaplan-Meier interval tables, multi-stratum
  contingency tables, reconstructed IPD, and linked migrant cohort
  datasets.
- [`plot_progression()`](https://n8thangreen.github.io/EvidSynthTB/reference/plot_progression.md)
  updated to support both uncured (`"uncured"`) and population mixture
  (`"population"`) survival curves from
  [`posterior::draws_df`](https://mc-stan.org/posterior/reference/draws_df.html)
  output.
