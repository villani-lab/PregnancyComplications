# Builds a Patient_ID -> Condition_group lookup table, the patient-level counterpart of
# "/projects/Pregnancy_placenta/Results/13_DE/sample_condition_list.csv" (which is keyed by the
# raw sequencing "Sample" ID instead). Source of truth: the one-hot Condition_* boolean obs
# columns baked into the GEO-release object by Misc_GEO_objects.py (lines ~152-161, 444-449),
# which were themselves derived from sample_condition_list.csv - so this file should agree with
# it, just re-keyed to Patient_ID and saved for reuse by any patient-level analysis.
# A patient can carry more than one condition flag (e.g. PrePE + FGR), so the output is long
# format: one row per (Patient_ID, Condition_group) membership, same many-to-many convention as
# sample_condition_list.csv (see Miscellaneous.py docstring).

import muon as mu
import pandas as pd

ALLCELLS_PATH = "/projects/Pregnancy_placenta/Results/21_Final_objects/GEO_objects/AllCells.h5mu"
OUT_PATH = "/projects/Pregnancy_placenta/Code/Publication_code/AuxFiles/patient_condition_list.csv"

# obs boolean column -> canonical Condition_group label (same 8 labels used everywhere else in
# the repo, e.g. Functions.py's condition_colors / 24_Differential_abundance.py comparisons)
CONDITION_FLAG_COLS = {
    "Condition_UnlaboredCont": "Cont_unlabored",
    "Condition_LaboredCont": "Cont_labored",
    "Condition_PTBCont": "Cont_PTB",
    "Condition_sPTB": "PTB",
    "Condition_PrePE": "PrePE",
    "Condition_PE": "TermPE",
    "Condition_FGR": "FGR",
    "Condition_DM1": "DM1",
}

data = mu.read(ALLCELLS_PATH)
obs = data["GEX"].obs
pt_flags = obs[["Patient_ID"] + list(CONDITION_FLAG_COLS.keys())].drop_duplicates().set_index("Patient_ID")

# Sanity check before trusting the mapping below - confirm these are really boolean-ish flags
print(pt_flags.dtypes)
print(pt_flags.apply(lambda c: c.unique()))

# Membership by exact match rather than bool()/astype(bool), since bool("False") is True in
# Python - this covers True/False, 1/0, or "x" flags regardless of which representation it turns
# out to be.
TRUTHY = {True, "True", "TRUE", 1, "1", "x", "X"}
flag_counts = pt_flags[list(CONDITION_FLAG_COLS.keys())].isin(TRUTHY).sum(axis=1)
print("Patients with zero condition flags set:", (flag_counts == 0).sum())
print("Patients with >1 condition flags set (e.g. PrePE+FGR):", (flag_counts > 1).sum())

rows = [(pt, label)
        for flag_col, label in CONDITION_FLAG_COLS.items()
        for pt in pt_flags.index[pt_flags[flag_col].isin(TRUTHY)]]
patient_condition = pd.DataFrame(rows, columns=["Patient_ID", "Condition_group"]).drop_duplicates()
patient_condition.to_csv(OUT_PATH, index=False)
print(patient_condition["Condition_group"].value_counts())
