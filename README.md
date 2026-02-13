# selex

Selectivity-method documentation for age-structured stock assessment models, authored in Quarto and rendered to standalone HTML files in `docs/`.

## What this repo contains

- `index.qmd`: main manuscript page that embeds component sections
- `logistic_selectivity.qmd`: standard logistic selectivity
- `double_logistic_selectivity.qmd`: 3-parameter double-logistic selectivity
- `spline_selectivity.qmd`: penalized spline selectivity
- `timevarying_selectivity.qmd`: time-varying selectivity with separable AR1 structure
- `R/`: reusable R helper functions
- `docs/`: rendered output for GitHub Pages

## Requirements

- Quarto CLI
- R
- R packages used in the notebooks:
  - `ggplot2`
  - `dplyr`
  - `tidyr`
  - `patchwork`
  - `splines` (base-recommended)
  - `RTMB`
  - `SparseNUTS` (for MCMC sections)

## Render standalone HTML

This project is configured in `_quarto.yml` to render the main QMD files to standalone HTML (`embed-resources: true`) into `docs/`.

```bash
quarto render
```

If you need to force re-execution of frozen chunks:

```bash
quarto render --cache-refresh
```

## Output files

Primary rendered pages:

- `docs/index.html`
- `docs/logistic_selectivity.html`
- `docs/double_logistic_selectivity.html`
- `docs/spline_selectivity.html`
- `docs/timevarying_selectivity.html`

These are self-contained HTML files suitable for direct sharing or static hosting.

## Local preview

```bash
quarto preview
```

## Publishing

The repository is set up to serve content from `docs/`. Typical publish flow:

1. Render with `quarto render`.
2. Commit source and updated `docs/` artifacts.
3. Push to `main`.

## Notes on reproducibility

- Execution freezing is enabled (`execute.freeze: true`) to keep renders stable.
- Use `--cache-refresh` when code/data changes require recomputation.
