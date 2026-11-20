# EvidSynthTB 0.2.0

* Migrated Bayesian backend to `cmdstanr` and `posterior` with modernized Stan models in `inst/stan/` (`ltbi_module.stan` and `evidsynth_modular.stan`).
* `compare_evidsynth()`, `plot_model_comparison()`, and `plot_survival_comparison()` added to fit and compare multiple synthesis strategies side-by-side.
* `create_fake_data()` and `create_fake_published_data()` implemented to simulate individual-level cohorts, published subgroup summary tables, and UKHSA-style annual surveillance progression tables with imperfect tests (`sens`, `spec`), covariates, and prior/data conflict scenarios.
* `diagnose_feedback()` and `evaluate_comparison()` added to quantify inter-module conflict, feedback shift relative to the cut model, survival WAIC, parameter bias/RMSE/coverage, and survival curve integrated absolute error.
* `evidsynth_fit()` expanded to natively accept both published aggregate/interval/KM tables and individual-level survival data across six modular synthesis strategies (`"joint"`, `"cut"`, `"independent"`, `"melding"`, `"elastic"`, `"tempered"`), `"population"` and `"latent_class"` coupling, `"all"` and `"pos_only"` progression cohorts, and `"gompertz"`, `"weibull"`, and `"exponential"` survival distributions.
* `load_abubakar_km()`, `load_predict_aggregated()`, `load_abubakar_2018()`, and `load_predict_migrant()` added to load published Kaplan-Meier interval tables, multi-stratum contingency tables, reconstructed IPD, and linked migrant cohort datasets.
* `plot_progression()` updated to support both uncured (`"uncured"`) and population mixture (`"population"`) survival curves from `posterior::draws_df` output.
