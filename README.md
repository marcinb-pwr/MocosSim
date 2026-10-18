# MocosSim

[![Build Status](https://app.travis-ci.com/MOCOS-COVID19/MocosSim.svg?branch=master)](https://app.travis-ci.com/MOCOS-COVID19/MocosSim)
[![Test Coverage](https://codecov.io/github/MOCOS-COVID19/MocosSim/coverage.svg?branch=master)](https://codecov.io/github/MOCOS-COVID19/MocosSim?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/MOCOS-COVID19/MocosSim/badge.svg?branch=master)](https://coveralls.io/github/MOCOS-COVID19/MocosSim?branch=master)

Simulation engine for coronavirus spreading

## Reinfections and Omicron cross-immunity

The simulator keeps every infection as a separate episode, so a recovered
person can be infected again during the same run. Cross-immunity is configured
when loading or making parameters. Rows of `infection_cross_immunity` identify
the previous strain and columns identify the challenge strain, in the order
Chinese, British, Delta, Omicron. Values are probabilities in `[0, 1]`.

```julia
params = MocosSim.load_params(
  population=population,
  # Defaults use the literature-calibrated tables. To run a sensitivity draw:
  cross_immunity_std=MocosSim.calibrated_cross_immunity_std(),
  vaccination_immunity_std=0.05,
  # other scenario parameters...
)
```

The default infection table is `[0.90 0.88 0.85 0.45; 0.88 0.90 0.86 0.45;
0.82 0.84 0.90 0.50; 0.40 0.42 0.45 0.60]`. It is a peak
infection-protection proxy, not protection against severe disease. Rows are the
previous strain and columns are the challenge strain.
Protection wanes exponentially from the latest infection with a default
500-day half-life. Vaccination protection uses the latest vaccination date, a
default 100-day half-life, and the recent-mRNA-booster proxy
`[0.90, 0.90, 0.95, 0.67]`; it also becomes zero at the scheduled
loss-of-immunity event. A scenario modelling a primary series, a particular
vaccine product, or a different outcome should override these defaults.

`cross_immunity_std` accepts either one standard deviation or a 4×4 table;
`vaccination_immunity_std` accepts one value or a four-element vector. A
truncated-normal draw is made once per parameter set using the supplied RNG,
which makes it suitable for reproducible sensitivity ensembles without adding
infection-by-infection noise. A zero standard deviation gives the central
calibration exactly. The four `vaccination_cross_immunity` entries allow
immune escape to be calibrated independently for each challenge strain (for
example, lower vaccine protection against Omicron). Vaccine and prior-infection
protection are combined as independent layers, so hybrid immunity is supported.
Saved results now contain one row per infection and include
`infection_subjects`, allowing repeated subject IDs across the two-year run.

### Calibration references

The central values deliberately synthesize epidemiological effectiveness rather
than converting neutralization titres directly into probabilities. They are a
pragmatic model calibration: the source studies mix protection against any and
symptomatic infection, and reverse-chronology cells (for example Omicron
followed by Delta) were rarely observed. Those cells have deliberately larger
suggested sensitivity deviations and should not be interpreted as direct point
estimates. The standard-deviation table is a sensitivity range, not a published
pooled standard error:

* COVID-19 Forecasting Team, *Past SARS-CoV-2 infection protection against
  re-infection: a systematic review and meta-analysis*, **The Lancet** 2023,
  DOI [10.1016/S0140-6736(22)02465-5](https://doi.org/10.1016/S0140-6736(22)02465-5).
  This supports high (>82%) pre-Omicron protection and substantially lower
  protection against Omicron BA.1.
* Altarawneh et al., *Effects of Previous Infection and Vaccination on
  Symptomatic Omicron Infections*, **NEJM** 2022, DOI
  [10.1056/NEJMoa2203965](https://doi.org/10.1056/NEJMoa2203965). This supports
  approximately 46% pre-Omicron infection protection against BA.2 and stronger
  hybrid immunity.
* Andrews et al., *Covid-19 Vaccine Effectiveness against the Omicron Variant*,
  **NEJM** 2022, DOI
  [10.1056/NEJMoa2119451](https://doi.org/10.1056/NEJMoa2119451). This anchors
  the recent-booster Omicron value and the faster vaccine waning assumption.
* Reynolds et al., *Immune boosting by B.1.1.529 (Omicron) depends on previous
  SARS-CoV-2 exposure*, **Science** 2022, DOI
  [10.1126/science.abq1841](https://doi.org/10.1126/science.abq1841), together
  with Rössler et al., *Omicron infection enhances Delta antibody immunity in
  vaccinated persons*, **Nature** 2022, DOI
  [10.1038/s41586-022-04830-x](https://doi.org/10.1038/s41586-022-04830-x),
  motivates treating the reverse Omicron-to-pre-Omicron cells as more uncertain.
