# Simulate published summary tables and UKHSA-style annual surveillance tables

Generates synthetic aggregate tables matching the two most common forms
of published TB literature and public health surveillance reports
(without individual participant data):

1.  \`summary_table\`: Subgroup-level published summary table reporting
    group sample size (\`n\`), positive LTBI test count (\`y\`),
    cumulative active TB case count (\`d\`) over follow-up window \`\[0,
    t_horizon\]\`, and subgroup covariates.

2.  \`interval_table\`: Annual surveillance report table (e.g., UKHSA /
    ETS reports) giving active TB case counts (\`d\`) in annual
    intervals \`\[t_start, t_end\]\` (Year 0-1, 1-2, ..., \`t_horizon\`)
    out of an initial screened / migrant cohort (\`n\`).

## Usage

``` r
create_fake_published_data(
  n_groups = 6,
  n_per_group = 500,
  n_surveillance = 4000,
  p_ltbi = 0.3,
  p_ltbi_surv = NULL,
  sens = 0.85,
  spec = 0.95,
  surv_dist = c("gompertz", "weibull", "exponential"),
  shape = -0.3,
  rate = 0.06,
  t_horizon = 5,
  seed = NULL
)
```

## Arguments

- n_groups:

  Integer. Number of published demographic/clinical strata in
  \`summary_table\`. Default is \`6\`.

- n_per_group:

  Integer. Sample size per group in \`summary_table\`. Default is
  \`500\`.

- n_surveillance:

  Integer. Total denominator cohort size for the annual surveillance
  \`interval_table\`. Default is \`4000\`.

- p_ltbi:

  Numeric in \`(0, 1)\`. True baseline LTBI prevalence in the published
  LTBI survey table. Default is \`0.30\`.

- p_ltbi_surv:

  Numeric in \`(0, 1)\` or \`NULL\`. True LTBI prevalence in the
  surveillance progression cohort. Defaults to \`p_ltbi\` (set different
  to simulate conflict / under-reporting / case-mix shift between
  published surveys and routine surveillance).

- sens:

  Numeric in \`(0, 1\]\`. Diagnostic test sensitivity. Default \`0.85\`.

- spec:

  Numeric in \`(0, 1\]\`. Diagnostic test specificity. Default \`0.95\`.

- surv_dist:

  Character. Uncured survival family: \`"gompertz"\`, \`"weibull"\`, or
  \`"exponential"\`.

- shape:

  Numeric. Shape parameter for Gompertz or Weibull. Default \`-0.3\`
  (decaying progression hazard over years since entry).

- rate:

  Numeric (\`\> 0\`). Initial annual progression hazard rate \`lambda\`.
  Default \`0.06\`.

- t_horizon:

  Integer (\`\>= 1\`). Number of annual follow-up intervals. Default is
  \`5\`.

- seed:

  Optional integer random seed.

## Value

A list with \`summary_table\` (cumulative subgroup counts),
\`interval_table\` (annual interval counts by time since entry),
\`prevalence_data\` (aggregate survey totals), and \`truth\` (true
parameters).
