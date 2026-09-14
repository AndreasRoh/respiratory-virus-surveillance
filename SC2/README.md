# SARS-CoV-2: normal analysis

Run this workflow from the repository root. The combined analysis driver controls the dependency order, so run only the command below rather than sourcing individual files first.

```r
source("SC2/SC2_Analysis.R")
```

Execution order inside the driver:

1. `SC2_SQLquery_BNCOVID19.R` — load the historical BNCOVID19 extract.
2. `SC2_DataCleaning_BNCOVID19.R` — clean the historical extract.
3. `SC2_SQLquery_25-26.R` — load the current-season extract.
4. `SC2_DataCleaning.R` — clean and combine the current-season data.
5. `SC2_Classification.R` — add Pangolin, Tessy, and origin classifications.
6. `SC2_Analysis.R` — create the normal surveillance report and exports.
