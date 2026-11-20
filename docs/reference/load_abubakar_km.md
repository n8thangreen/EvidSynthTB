# Load published Kaplan-Meier interval life table from Abubakar et al. (2018)

Loads the digitized Kaplan-Meier interval table (\`KMdata.txt\`,
containing 69 time intervals with \`n_risk\` and \`d\` progression
events) together with the published aggregate baseline IGRA positivity
count (\`y = 1444\` positive out of \`n = 6380\` tested). Fitting
\[evidsynth_fit()\] to this published life table runs orders of
magnitude faster than fitting to 6,380 reconstructed individual rows.

## Usage

``` r
load_abubakar_km()
```

## Value

A list with \`progression_data\` (tibble with \`t_start\`, \`t_end\`,
\`n_risk\`, \`d\`, \`n_censored\`) and \`prevalence_data\` (tibble with
\`y = 1444L\`, \`n = 6380L\`, \`test_id = 1L\`, \`mu_hat\`,
\`sigma_hat\`).
