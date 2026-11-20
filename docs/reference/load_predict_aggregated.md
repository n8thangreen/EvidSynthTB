# Load published-style aggregated PREDICT-TB subgroup contingency table

Loads the grouped contingency table (\`aggregated_data.csv\`)
summarizing sample size (\`pop\`), positive LTBI test counts (\`ltbi\`),
and cumulative active TB progression counts (\`tb\`) stratified by
\`age_grp\`, \`sex\`, \`ethnicity\`, TB incidence in country of birth
(\`inc_cob_participant2\`), years since UK entry
(\`yearssinceentry_grp\`), prior BCG vaccination (\`prevbcg\`), and
\`reasonforscreening\` (\`"Contact"\` vs \`"Migrant"\`).

## Usage

``` r
load_predict_aggregated(
  collapse_by = c("reasonforscreening", "inc_cob_participant2", "age_grp"),
  t_horizon = 5
)
```

## Arguments

- collapse_by:

  Optional character vector of grouping columns to collapse over
  (summing \`pop\`, \`ltbi\`, and \`tb\`). Defaults to
  \`c("reasonforscreening", "inc_cob_participant2", "age_grp")\`. Set to
  \`NULL\` to return all raw cells.

- t_horizon:

  Numeric (\`\> 0\`). Cumulative follow-up horizon in years assigned to
  \`t_end\`. Default is \`5\`.

## Value

A list with \`summary_table\` (tibble with \`n\`, \`y\`, \`d\`,
\`t_start\`, \`t_end\`, and grouping covariates) and \`prevalence_data\`
(same table ready for \[evidsynth_fit()\]).
