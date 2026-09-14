# RSV: normal analysis

Run this workflow from the repository root. The analysis driver controls the dependency order, so run only the command below rather than sourcing individual files first.

```r
source("RSV/RSV_Analysis.R")
```

Execution order inside the driver:

1. `RSV_SQLquery.R` — load the RSV extract.
2. `RSV_DataCleaning.R` — clean and prepare `rsvdb`.
3. `RSV_Analysis.R` — create the normal surveillance report and exports.
