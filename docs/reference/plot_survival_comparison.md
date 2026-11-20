# Plot survival curve comparisons across synthesis models

Plots posterior median and 95 curve across multiple fitted synthesis
strategies from \[compare_evidsynth()\].

## Usage

``` r
plot_survival_comparison(comparison_obj, type = c("uncured", "population"))
```

## Arguments

- comparison_obj:

  An \`"evidsynth_comparison"\` object returned by
  \[compare_evidsynth()\].

- type:

  Character. \`"uncured"\` (\`S_pred\`) or \`"population"\`
  (\`S_pop_pred\`).

## Value

A \`ggplot2\` object.
