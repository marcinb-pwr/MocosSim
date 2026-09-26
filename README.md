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
cross_immunity = [
  0.90 0.80 0.70 0.15
  0.80 0.90 0.75 0.20
  0.70 0.75 0.90 0.25
  0.20 0.25 0.30 0.80
]

params = load_params(
  population=population,
  infection_cross_immunity=cross_immunity,
  vaccination_cross_immunity=[0.90, 0.85, 0.80, 0.35],
  cross_immunity_half_life=180.0,
  vaccination_immunity_half_life=120.0,
  # other scenario parameters...
)
```

These numbers are an example, not a calibration. Scenario-specific estimates
should be supplied explicitly. Protection wanes exponentially from the time of
the latest infection; use `cross_immunity_half_life=Inf` to disable that
waning. Vaccination protection uses the time of the latest vaccination and its
own `vaccination_immunity_half_life`; it also becomes zero at the scheduled
loss-of-immunity event. The four `vaccination_cross_immunity` entries allow
immune escape to be calibrated independently for each challenge strain (for
example, lower vaccine protection against Omicron). Vaccine and prior-infection
protection are combined as independent layers, so hybrid immunity is supported.
Saved results now contain one row per infection and include
`infection_subjects`, allowing repeated subject IDs across the two-year run.
