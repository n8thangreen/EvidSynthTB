# Load cleaned individual-level PREDICT-TB migrant cohort data

Loads the individual-level migrant sub-cohort from
\`cleaned_migrant_predict_data.RData\`, containing baseline diagnostic
test positivity (\`x\`), follow-up time (\`t\`), active TB progression
indicator (\`d\`), centered age (\`age_z\`), and \`ethnicity\`.

## Usage

``` r
load_predict_migrant()
```

## Value

A list with \`progression_data\` (tibble with individual \`t\`, \`d\`,
\`x\`, \`age\`, \`age_z\`, \`ethnicity\`) and \`prevalence_data\`
(aggregate summary of test positivity in the cohort).
