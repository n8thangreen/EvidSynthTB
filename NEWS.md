# EvidSynthTB 0.2.0

* Migrated Bayesian backend to `cmdstanr` and `posterior` with modernized Stan models in `inst/stan/` (`ltbi_module.stan` and `evidsynth_modular.stan`).
* `compare_evidsynth()`, `plot_model_comparison()`, and `plot_survival_comparison()` added to fit and compare multiple synthesis strategies side-by-side.
* `create_fake_data()` implemented to simulate linked or unlinked LTBI diagnostic and active TB progression survival cohorts with imperfect tests (`sens`, `spec`), multiple tests per person, `pos_only` subgroup follow-up, covariates, and prior/data conflict scenarios.
* `evidsynth_fit()` expanded to support six modular synthesis strategies (`"joint"`, `"cut"`, `"independent"`, `"melding"`, `"elastic"`, `"tempered"`), both `"population"` and `"latent_class"` coupling, `"all"` and `"pos_only"` progression cohorts, and `"gompertz"`, `"weibull"`, and `"exponential"` survival distributions.
* `load_abubakar_2018()` and `load_predict_migrant()` added to load the reconstructed UK PREDICT-TB Kaplan-Meier dataset and cleaned migrant cohort dataset.
* `plot_progression()` updated to support both uncured (`"uncured"`) and population mixture (`"population"`) survival curves from `posterior::draws_df` output.
