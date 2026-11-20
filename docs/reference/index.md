# Package index

## Model Fitting & Multi-Model Synthesis

Fit and compare Bayesian Multi-Parameter Evidence Synthesis (MPES)
models in CmdStan.

- [`evidsynth_fit()`](https://n8thangreen.github.io/EvidSynthTB/reference/evidsynth_fit.md)
  : Fit Bayesian Multi-Parameter Evidence Synthesis Model for LTBI and
  TB Progression
- [`compare_evidsynth()`](https://n8thangreen.github.io/EvidSynthTB/reference/compare_evidsynth.md)
  : Compare multiple Bayesian evidence synthesis models

## Quantitative Diagnostics & Visualization

Evaluate parameter recovery, diagnose inter-module conflict and feedback
shift, and plot posterior comparisons.

- [`evaluate_comparison()`](https://n8thangreen.github.io/EvidSynthTB/reference/evaluate_comparison.md)
  : Evaluate parameter recovery and survival curve error against ground
  truth
- [`diagnose_feedback()`](https://n8thangreen.github.io/EvidSynthTB/reference/diagnose_feedback.md)
  : Diagnose inter-module conflict, feedback shift, and survival WAIC
- [`plot_progression()`](https://n8thangreen.github.io/EvidSynthTB/reference/plot_progression.md)
  : Plot posterior progression survival curve
- [`plot_model_comparison()`](https://n8thangreen.github.io/EvidSynthTB/reference/plot_model_comparison.md)
  : Plot posterior parameter comparisons across synthesis models
- [`plot_survival_comparison()`](https://n8thangreen.github.io/EvidSynthTB/reference/plot_survival_comparison.md)
  : Plot survival curve comparisons across synthesis models

## Data Simulation & Real TB Datasets

Simulate synthetic LTBI/TB datasets (IPD or published tables) and load
UK PREDICT-TB datasets.

- [`create_fake_data()`](https://n8thangreen.github.io/EvidSynthTB/reference/create_fake_data.md)
  : Simulate synthetic LTBI diagnostic and TB progression survival data
- [`create_fake_published_data()`](https://n8thangreen.github.io/EvidSynthTB/reference/create_fake_published_data.md)
  : Simulate published summary tables and UKHSA-style annual
  surveillance tables
- [`load_abubakar_km()`](https://n8thangreen.github.io/EvidSynthTB/reference/load_abubakar_km.md)
  : Load published Kaplan-Meier interval life table from Abubakar et al.
  (2018)
- [`load_predict_aggregated()`](https://n8thangreen.github.io/EvidSynthTB/reference/load_predict_aggregated.md)
  : Load published-style aggregated PREDICT-TB subgroup contingency
  table
- [`load_abubakar_2018()`](https://n8thangreen.github.io/EvidSynthTB/reference/load_abubakar_2018.md)
  : Load reconstructed Abubakar (2018) UK PREDICT-TB data
- [`load_predict_migrant()`](https://n8thangreen.github.io/EvidSynthTB/reference/load_predict_migrant.md)
  : Load cleaned individual-level PREDICT-TB migrant cohort data
