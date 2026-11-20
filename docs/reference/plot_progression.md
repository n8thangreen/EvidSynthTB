# Plot posterior progression survival curve

Plots the posterior median, mean, 50 interval of the survival curve over
time from an \`evidsynth_fit\` object.

## Usage

``` r
plot_progression(evidsynth_obj, type = c("uncured", "population"))
```

## Arguments

- evidsynth_obj:

  Output of \[evidsynth_fit()\].

- type:

  Character. Either \`"uncured"\` (uncured progression survival
  \`S_pred\`, conditional on being LTBI positive) or \`"population"\`
  (overall mixture survival \`S_pop_pred\`, including the cured fraction
  \`1 - prev_cf\`).

## Value

A \`ggplot2\` object.
