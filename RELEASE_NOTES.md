# Release notes

## v1.1.0

Production release prepared from the verified develop baseline.

### Repository layout

- main is the production branch.
- develop is the development branch.
- INF, RSV, and SC2 contain their executable modules.
- Source_files contains shared runtime helpers and reference inputs.
- Generated output is intentionally excluded from Git.

### Naming convention

- Executable pathogen modules use PATHOGEN_Role.R.
- Analysis, DataCleaning, Statistics, QualityControl, and export roles use title case.
- Database-specific SQL query and BNCOVID19 data-cleaning filenames are retained.
- Shared helpers use lower snake case in Source_files.

### Branch separation

The former in-repository development snapshot was removed. Each Git branch now
contains its own complete codebase and has no cross-branch runtime dependency.
