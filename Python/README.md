# Python code

Code that generates the figures of the paper that were made in Python. R code is in `../R code`.

## Layout

| Folder / file | Content |
|---|---|
| `Figure_*.py`, `Figure_*.ipynb` | Code for main figure panels (named by figure/panel) |
| `SuppFig_*.py`, `SuppFig_*.ipynb` | Code for supplementary figures |
| `SuppTable_*.ipynb` | Code for supplementary tables |
| `Functions/` | Helper functions imported by the figure code |
| `AuxFiles/` | Auxiliary input files (marker gene lists, condition lists, thresholds, ...) |
| `environments/` | Package lists of the environments used (see below) |

## Setup: set `code_dir`

Every script/notebook starts with

```python
code_dir="/projects/Pregnancy_placenta/Code_Publication/Python/"
```

Change this single variable to the location of this folder (keep the trailing `/`). Data and output paths are not yet parameterised.

## Environments

The analyses were run in different environments:

| Environment | Used for | File |
|---|---|---|
| CellTypist | SuppFig_2 (CellTypist comparison) | `environments/celltypist_env.txt` |
| pregnancy_placenta (general single-cell; Python 3.8) | Figure_1/3/4/5/6 scripts and notebooks, SuppFig_1/7/8/9/10 | `environments/pregnancy_placenta_env.txt` |
| spatial_env (Visium) | Figure_1H-I, Figure_2A-F, Figure_6 + SuppFig_10 Visium notebook, SuppFig_4, SuppFig_5, SuppTable_6 | *to be added* |
| SuppFig_11 (Miltenyi MACSima / Cellpose) | SuppFig_11.py | *to be added* |

Both files are `pip freeze` lists: `pip install -r environments/<file>`. `pregnancy_placenta_env.txt` includes the `pseudobulk` package installed from GitHub (villani-lab/pseudobulk, pinned commit) and GPU builds of torch (CUDA 12.1).

> The env/file assignments above are provisional; they will be confirmed when the environment files are added.
