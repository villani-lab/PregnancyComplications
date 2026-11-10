import os
import re
import glob
import pandas as pd

def parse_ipa_filename(path):
    """
    Expected: {cell_type}-{condition}_vs_{control}-{annotation}.xls(x)
    Example: HBC_11-PTB_vs_Cont_unlabored-diseases.xls
    """
    base = os.path.basename(path)

    pattern = (
        r"(?P<cell_type>.+?)"
        r"-(?P<condition>.+?)_vs_(?P<control>.+?)"
        r"-(?P<annotation>.+?)\.(?P<ext>xls|xlsx|csv|tsv)$"
    )

    m = re.match(pattern, base, flags=re.IGNORECASE)
    if not m:
        raise ValueError(f"Filename does not match expected pattern: {base}")

    d = m.groupdict()
    return {
        "cell_type": d["cell_type"],
        "condition": d["condition"],
        "control": d["control"],
        "annotation": d["annotation"],
        "ext": d["ext"].lower()
    }

def read_one_ipa_file(path, skiprows=1):
    meta = parse_ipa_filename(path)
    ext = meta["ext"]

    if ext in ("xls", "xlsx"):
        df = pd.read_excel(path, skiprows=skiprows)
    elif ext == "csv":
        df = pd.read_csv(path)
    elif ext == "tsv":
        df = pd.read_csv(path, sep="\t")
    else:
        raise ValueError(f"Unsupported extension: {ext}")

    # Add metadata columns
    df["cell_type"] = meta["cell_type"]
    df["condition"] = meta["condition"]
    df["control"] = meta["control"]
    df["annotation_type"] = meta["annotation"]

    return df

def read_and_merge_IPA_outputs(loc_of_IPA_gene_files, recursive=False, skiprows=1, verbose=True):
    patterns = ["*.xls", "*.xlsx", "*.csv", "*.tsv"]
    all_files = []
    for pat in patterns:
        if recursive:
            all_files.extend(glob.glob(os.path.join(loc_of_IPA_gene_files, "**", pat), recursive=True))
        else:
            all_files.extend(glob.glob(os.path.join(loc_of_IPA_gene_files, pat)))

    all_files = sorted(set(all_files))

    if verbose:
        print(f"Found {len(all_files)} files in: {loc_of_IPA_gene_files}")
        if len(all_files) > 0:
            print("Example files:", [os.path.basename(x) for x in all_files[:5]])

    dfs = []
    failed = []

    for f in all_files:
        try:
            df = read_one_ipa_file(f, skiprows=skiprows)
            dfs.append(df)
        except Exception as e:
            failed.append((f, str(e)))

    if verbose and failed:
        print(f"\nSkipped {len(failed)} files due to errors. First 5:")
        for f, err in failed[:5]:
            print(" -", os.path.basename(f), "->", err)

    if len(dfs) == 0:
        raise ValueError(
            "No objects to concatenate. "
            "Either no files were found, or all files failed parsing/reading.\n"
            f"Directory: {loc_of_IPA_gene_files}\n"
            f"Tip: set recursive=True if files are in subfolders."
        )

    return pd.concat(dfs, ignore_index=True)
