# Influenza: normal analysis

Run this workflow from the repository root. The analysis driver controls the dependency order, so run only the command below rather than sourcing individual files first.

```r
source("INF/INF_Analysis.R")
```

Execution order inside the driver:

1. `INF_SQLquery_25-26.R` — load the current influenza extract.
2. `INF_DataCleaning.R` — clean and prepare `fludb`.
3. `INF_QualityControl.R` — run data-quality checks.
4. `INF_Analysis.R` — create the normal surveillance report and exports.

Excluded workflows: `INF_TESSyExport.R` and `INF_RAVNExport.R`.
