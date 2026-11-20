# Load reconstructed Abubakar (2018) UK PREDICT-TB data

Loads the digitized individual participant survival data
(\`IPDdata.txt\`) reconstructed from the Kaplan-Meier curve in Abubakar
et al. (2018) alongside the aggregate baseline IGRA positivity count
(\`y = 1444\` positive out of \`n = 6380\` tested).

## Usage

``` r
load_abubakar_2018()
```

## Value

A list with \`progression_data\` (tibble with columns \`t\` and \`d\`)
and \`prevalence_data\` (tibble with \`y = 1444L\`, \`n = 6380L\`,
\`test_id = 1L\`, and empirical logit plug-in estimates \`mu_hat\` and
\`sigma_hat\`).
