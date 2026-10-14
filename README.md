# Sagittal fractal-dimension survival analysis

Code for extracting two- and three-dimensional box-counting fractal-dimension features from NIfTI tumor masks and evaluating survival models with a training set and two validation sets.

## Repository contents

- `src/fractal_dimension.py`: reusable command-line feature extractor.
- `notebooks/fractal_dimension_extraction.ipynb`: interactive feature-extraction example.
- `analysis/survival_analysis.R`: survival modeling, validation, calibration, and explainability analyses.
- `R/evaluation_helpers.R`: shared survival-model evaluation functions.

Clinical records, imaging masks, derived patient-level tables, and model outputs are not included because they may contain sensitive data.

## Python setup and feature extraction

Python 3.10 or newer is recommended.

```bash
python -m venv .venv
python -m pip install -r requirements.txt
python src/fractal_dimension.py \
  --input-dir data/masks \
  --output-csv results/fractal_dimension_features.csv
```

The extractor converts each NIfTI image to the closest canonical orientation before calculating axial, coronal, sagittal, and whole-volume features. Use `--seed`, `--repeats`, and `--shifts` to control reproducibility and computation time.

## R survival analysis

The R script reads paths from environment variables. If variables are omitted, it expects the following files under `data/`:

| Environment variable | Default file |
| --- | --- |
| `SAGITTAL_FD_DEVELOPMENT_CLINICAL` | `data/development_clinical.sav` |
| `SAGITTAL_FD_DEVELOPMENT_FEATURES` | `data/development_fd_features.csv` |
| `SAGITTAL_FD_EXTERNAL_CLINICAL` | `data/external_clinical.xlsx` |
| `SAGITTAL_FD_EXTERNAL_FEATURES` | `data/external_fd_features.csv` |
| `SAGITTAL_FD_OUTPUT_DIR` | `results/` |

Run the analysis from the repository root:

```bash
Rscript analysis/survival_analysis.R
```

The script checks for required input files and R packages before starting. Generated figures, fitted models, intermediate objects, and patient-level tables are written to `results/`, which is excluded from version control.

## Data naming

The analysis uses `train_data`, `validation_data1`, and `validation_data2`. Figures and result tables use the corresponding labels “Validation Set 1” and “Validation Set 2”. The `external_*` input-path variables are retained to document that the second source cohort was collected externally.

## License

Released under the MIT License. See `LICENSE`.
