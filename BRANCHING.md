# Branch Workflow

This repository uses two long-lived branches:

- `main`: production/stable branch
- `develop`: integration and testing branch

## Daily Workflow

- Start new feature or fix branches from `develop`.
- Merge completed work back into `develop` after testing.
- Avoid direct development on `main`.

## Promotion to Production

- Freeze the release candidate on `develop`.
- Run the relevant smoke tests:
  - `INF/INF_Analysis.R` when INF-related code changed
  - `RSV/RSV_Analysis.R` when RSV-related code changed
  - `SC2/SC2_Analysis_combined_singlefile.R` when SC2-related code changed
- Verify output artifacts open correctly for the changed pathogen(s), especially PPT export.
- Merge the exact tested `develop` commit into `main` with no additional code changes.

## Notes

- Shared helper changes in `Source_files` require smoke checks for all three pathogens before promotion.
- `main` should be treated as the production baseline.
- `develop` is the default branch for active integration work.
