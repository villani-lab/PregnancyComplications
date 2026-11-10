# Python code reorganisation – change log

Only file names/locations and comment/markdown header text were changed. No executable code was edited.

## File renames (old figure number -> revised-paper number)
| Old | New | Reason |
|---|---|---|
| SuppFig_7.py | SuppFig_10.py | TCR and tolerance (S7 -> S10) |
| SuppFig_6.py | SuppFig_9.py | DM1 and FGR (S6 -> S9) |
| SuppFig_5.py | SuppFig_7.py | Cellular dynamics vs labored controls (S5 -> S7) |
| SuppFig_4.ipynb | SuppFig_5.ipynb | Marker genes on visium cores (S4 -> S5) |
| SuppFig_3.ipynb | SuppFig_4.ipynb | Overview of visium cores (S3 -> S4) |
| Visium_publication_code/SuppFig_3.ipynb | Visium_publication_code/SuppFig_4.ipynb | byte-identical duplicate of the above (see open questions) |
| Figure_6_and_SuppFig_7-visium_plots.ipynb | Figure_6_and_SuppFig_10-visium_plots.ipynb | contains Supp 7K-M -> 10K-M |

Unchanged: Figure_1..6 files, SuppFig_1.py (S1 stays S1).

## Header (comment / markdown) edits
- SuppFig_7.py: `supp fig 5` -> `supp fig 7`; `Sup Fig 5A)`, `5B)` -> `7A)`, `7B)`
- SuppFig_9.py: `SupFig 6 A-C` -> `SupFig 9 A-C`
- SuppFig_10.py: `Sup Fig 7A)...7I, J)` -> `Sup Fig 10A)...10I, J)` (9 headers)
- SuppFig_5.ipynb: `Sup Fig4` -> `Sup Fig5`; `Sup Fig 4A-D` -> `Sup Fig 5A-D`
- SuppFig_4.ipynb (+ the Visium_publication_code copy): `SupFig3 - Visium` -> `SupFig4 - Visium`
- Figure_6_and_SuppFig_10-visium_plots.ipynb: `Sup Fig 7K`, `7L/M`, `7L`, `7M` -> `10K`, `10L/M`, `10L`, `10M`

## New figures added from "Code to add" (Python parts)
| From | To | Header |
|---|---|---|
| S2_celltypist_clean.py | SuppFig_2_celltypist_clean.py | added 1-line header comment at top |
| S8_myeloid_patient_level.py | SuppFig_8_myeloid_patient_level.py | added 1-line header comment at top |
| figS11_code/ (4 files) | SuppFig_11_code/ | folder renamed only; file names kept (see below) |

SuppFig_11_code keeps `run_figs11.py`, `figs11_*_helpers.py` and `thresholds.json` unchanged because the code imports them by name / loads the JSON beside the script, and `figs11_segmentation_helpers.py` references `run_figs11.py` by name.

## Deliberately NOT changed (code strings, not headers)
Output file names and paths inside the code still carry old numbers, e.g. `SuppFig6B_...pdf` in SuppFig_10.py, `Sup_Fig3A_...pdf` in SuppFig_4.ipynb, `.../SuppFigure_5/` work dirs.

## Still pending / open
- R parts of new S2 (S2_celltypist_plot_dotplot.R) and S6 (figS6_*.R) remain in "Code to add" for the R step.
- Table-numbering references (SuppTable_4.ipynb, "Sup Table 11/12" in Figure_4D-F.ipynb, "Table S5" comment in Figure_3.py) not touched.

## Table renumbering (applied afterwards)
Code labels were off by one vs. the first-submission paper numbering, so I mapped by content: old paper number -> new paper number.
- SuppTable_4.ipynb -> SuppTable_6.ipynb (Visium HD QC metrics: old S5 -> new S6)
- Figure_4D-F.ipynb markdown: `Sup Table 11 - angiogenesis score` -> `Sup Table 13` (old S12 -> new S13)
- Figure_3.py comment on `final.to_excel("Table_S6.xlsx")`: `# Table S5` -> `# Table S7` (cell frequency and DA: old S6 -> new S7)
Not changed (code strings): `Table_S6.xlsx`, `Supp_Table_11_AngioGenesis_Score.csv`, `Supp_Table_12_mixed_model_analysis.csv`.
Not numbered in the code, so no header to fix: the mixed-model cell ("Sup Table - Angio Score + Mixed MOdel", new S14), and the gene-table exports in Figure_4A_C (new S11) and Figure_5A_C-E (new S15).

## Table output file names (code strings, changed at user's request)
- Figure_3.py: `Table_S6.xlsx` -> `Table_S7.xlsx`
- Figure_4D-F.ipynb: `Supp_Table_11_AngioGenesis_Score.csv` -> `Supp_Table_13_AngioGenesis_Score.csv`
- Figure_4D-F.ipynb: `Supp_Table_12_mixed_model_analysis.csv` -> `Supp_Table_14_mixed_model_analysis.csv`
(This supersedes the "Not changed (code strings)" note in the table section above for these three names.)
