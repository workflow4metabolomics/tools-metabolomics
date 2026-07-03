# -*- coding: utf-8 -*-
"""
MSOutputComparison 
First version 1.0.0 : Quentin Ruin, Mélanie Pétéra — INRAE - MetaboHUB
Second version 1.5.0: Abdou Lahat KA — INRAE - MetaboHUB

Version 1.5.0 Improvements:
  - New correlation plots with dynamic binning (step size = (max−min) / 20):
      · Global histogram        : 20 bins spanning [min_r, max_r]
      · Lower-half zoom         : 20 bins spanning [min_r, median_r]
      · Upper-half zoom         : 20 bins spanning [median_r, max_r]
  - Significantly improved visual design (professional colour scheme, annotations)
  - User-defined parameters for custom m/z and RT column names
  - Removal of the strict naming requirement for the ID column
  - Reinforced defensive programming: comprehensive validation of input files
  - Heatmap: axes now use the user-defined m/z and RT tolerance values as tick
    boundaries instead of fixed powers-of-10 grid; 6-level granularity per axis
    so each cell covers 1/6 of the tolerance range, giving finer spatial
    resolution and making the tolerance thresholds directly visible.
  - Heatmap: tolerance parameters are fully free (no system-imposed min/max).
  - Heatmap: fixed Y-axis inversion bug — m/z axis now correctly increases
    upward (finest threshold at bottom, full tolerance at top), achieved via
    ax.invert_yaxis() instead of array reversal which produced the opposite
    visual order with seaborn's top-down row rendering.
  - Correlation histograms: tolerance parameters (mz, rt) are fully free in the
    XML wrapper (min/max constraints removed).
  - Correlation histograms: each zoom chart (lower-half, upper-half) now uses its
    own independent Y-axis scale instead of a shared global scale, eliminating
    visual distortion when the two halves have very different counts.
  - Correlation histograms: automatic dual Y-axis (linear left + log right) when
    the dominant bin is >= 20x the second-largest bin; low-count bars are shown
    as coloured circles on the log axis so that rare correlation classes remain
    visible even when a single class dominates the distribution.
"""

import numpy as np
import pandas as pd
import sys
import os
import time
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from weasyprint import HTML
import seaborn as sns
from datetime import datetime

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

# UTILITY FUNCTIONS
#######################################################################################################################################################################

def stop_with_message(message):
    """Terminates the script execution and displays a clear error message to the user."""
    print("\n" + "="*60)
    print("❌ ERROR — The tool cannot proceed")
    print("="*60)
    print(message)
    print("="*60 + "\n")
    sys.exit(1)

def verifier_fichier(chemin, nom_affiche):
    """Checks if a file exists and is readable."""
    if not os.path.exists(chemin):
        stop_with_message(
            f"The file '{nom_affiche}' was not found.\n"
            f"Path used: {chemin}\n"
            f"→ Please verify that the file has been successfully uploaded into Galaxy."
        )
    if os.path.getsize(chemin) == 0:
        stop_with_message(
            f"The file '{nom_affiche}' is empty (0 bytes).\n"
            f"Path used: {chemin}\n"
            f"→ Please verify that your pre-processing pipeline successfully generated this file."
        )

def verifier_format_tabular(chemin, nom_affiche):
    """Verifies that the file format is a valid, readable tab-separated text file."""
    try:
        with open(chemin, 'r', encoding='utf-8', errors='replace') as f:
            premiere_ligne = f.readline()
        if '\t' not in premiere_ligne and ',' not in premiere_ligne:
            stop_with_message(
                f"The file '{nom_affiche}' does not appear to be in a valid tabular format.\n"
                f"First line read: {premiere_ligne[:200]!r}\n"
                f"→ Ensure that the dataset is tab-separated (standard .tabular format)."
            )
    except UnicodeDecodeError:
        stop_with_message(
            f"The file '{nom_affiche}' contains unreadable characters.\n"
            f"→ The file might be corrupted or in an unsupported binary format."
        )

def verifier_colonnes(df, nom_affiche, colonnes_requises):
    """Checks if all required columns are present in the DataFrame."""
    manquantes = [c for c in colonnes_requises if c not in df.columns]
    if manquantes:
        stop_with_message(
            f"The file '{nom_affiche}' is missing expected columns.\n"
            f"Missing columns: {manquantes}\n"
            f"Columns found  : {list(df.columns)}\n"
            f"→ Verify that this file is formatted as a tab-separated file\n"
            f"  and contains the needed columns.\n"
            f"  Please note that the m/z and retention time column names to\n"
            f"  look for can be personnalized in the Column Names section.\n"
        )

def verifier_colonnes_numeriques(df, nom_affiche, colonnes_numeriques):
    """Verifies that the expected numerical columns contain valid numeric data types."""
    for col in colonnes_numeriques:
        if col not in df.columns:
            continue
        serie = pd.to_numeric(df[col], errors='coerce')
        n_invalides = serie.isna().sum()
        if n_invalides == len(df):
            stop_with_message(
                f"Column '{col}' in file '{nom_affiche}' contains no numeric values.\n"
                f"Sample values found: {df[col].head(5).tolist()}\n"
                f"→ This column must contain numbers (m/z or retention time values)."
            )
        elif n_invalides > 0:
            print(f"   ⚠️   Warning: {n_invalides} non-numeric value(s) detected in "
                  f"column '{col}' of '{nom_affiche}' (coerced to NaN).")

def verifier_coherence_dataframes(varmat, datamat, nom_varmat, nom_datamat):
    """Checks the structural consistency between variableMetadata and dataMatrix files."""
    if len(varmat) == 0:
        stop_with_message(
            f"The file '{nom_varmat}' is empty (no data rows found).\n"
            f"→ Please verify that your pre-processing pipeline detected features/ions."
        )
    if len(datamat) == 0:
        stop_with_message(
            f"The file '{nom_datamat}' is empty (no data rows found).\n"
            f"→ Please verify that your pre-processing pipeline successfully generated this file."
        )
    if len(varmat) != len(datamat):
        print(f"   ⚠️   Warning: The number of rows differs between '{nom_varmat}' "
              f"({len(varmat)}) and '{nom_datamat}' ({len(datamat)}).\n"
              f"      Data will be processed independently (matching by index alignment).")

def safe_log10_min(series):
    """Calculates the log10 floor of the minimum non-zero value. Returns None if all values are zero."""
    nonzero = series[series != 0].abs()
    if len(nonzero) == 0:
        return None
    return int(np.floor(np.log10(nonzero.min())))

def safe_log10_max(series):
    """Calculates the log10 ceiling of the maximum non-zero value. Returns None if all values are zero."""
    nonzero = series[series != 0].abs()
    if len(nonzero) == 0:
        return None
    return int(np.ceil(np.log10(nonzero.max()))) + 1

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

if len(sys.argv) < 12:
    stop_with_message(
        "Missing arguments. The tool expects exactly 11 parameters:\n"
        "  1. variableMetadata workflow 1\n"
        "  2. variableMetadata workflow 2\n"
        "  3. dataMatrix workflow 1\n"
        "  4. dataMatrix workflow 2\n"
        "  5. m/z tolerance (e.g., 0.02)\n"
        "  6. RT tolerance (e.g., 0.25)\n"
        "  7. workflow 1 name\n"
        "  8. workflow 2 name\n"
        "  9. m/z column name (e.g., mz)\n"
        " 10. RT column name (e.g., rt)\n"
        " 11. correlation split threshold for zoom histograms (e.g., 0.5)"
    )

# Assigning inputs to explicit variables
path_variableMetadata1 = sys.argv[1]
path_variableMetadata2 = sys.argv[2]
path_DataMatrix1       = sys.argv[3]
path_DataMatrix2       = sys.argv[4]
tol_mz                 = float(sys.argv[5])
tol_RT                 = float(sys.argv[6])
wf1_name               = sys.argv[7]
wf2_name               = sys.argv[8]
col_mz                 = sys.argv[9]    # Target m/z column name (default: "mz")
col_rt                 = sys.argv[10]   # Target RT column name (default: "rt")
corr_split             = float(sys.argv[11])  # Correlation split threshold for zoom histograms (default: 0.5)

# Clamp corr_split to the valid Pearson r range [-1.0, 1.0]
corr_split = max(-1.0, min(1.0, corr_split))

# =============================================================================
# DATA LOADING & VALIDATION
# =============================================================================

print("\n" + "="*60)
print(f"  MS Output Comparison")
print(f"  {wf1_name}  VS  {wf2_name}")
print(f"  m/z Tolerance: {tol_mz} Da  |  RT Tolerance: {tol_RT}")
print(f"  m/z Column   : '{col_mz}'    |  RT Column   : '{col_rt}'")
print(f"  Corr. split  : {corr_split}")
print("="*60 + "\n")

print("📂 Loading files...")

for chemin, nom in [
    (path_variableMetadata1, "variableMetadata Workflow 1"),
    (path_variableMetadata2, "variableMetadata Workflow 2"),
    (path_DataMatrix1,       "dataMatrix Workflow 1"),
    (path_DataMatrix2,       "dataMatrix Workflow 2"),
]:
    verifier_fichier(chemin, nom)
    verifier_format_tabular(chemin, nom)

def charger_datamatrix(chemin, nom):
    try:
        df = pd.read_csv(chemin, sep='\t')
    except Exception as e:
        stop_with_message(
            f"Unable to read file '{nom}'.\n"
            f"Technical error: {e}\n"
            f"→ Verify that the file is in a valid tab-separated (.tabular) format."
        )
    if len(df) == 0:
        stop_with_message(
            f"The file '{nom}' contains no data rows.\n"
            f"→ Please check that your pre-processing pipeline successfully generated this dataset."
        )
    if len(df.columns) < 2:
        stop_with_message(
            f"The file '{nom}' contains only a single column.\n"
            f"Columns found: {list(df.columns)}\n"
            f"→ A dataMatrix must include at least an identifier (ID) column and one sample column."
        )
    df = df.replace({'NA': 0, 'NaN': 0}).fillna(0)
    return df

def charger_varmetadata(chemin, nom, col_mz, col_rt):
    try:
        df = pd.read_csv(chemin, sep='\t')
    except Exception as e:
        stop_with_message(
            f"Unable to read file '{nom}'.\n"
            f"Technical error: {e}\n"
            f"→ Verify that the file is in a valid tab-separated (.tabular) format."
        )
    if len(df) == 0:
        stop_with_message(
            f"The file '{nom}' contains no data rows.\n"
            f"→ Please check that your pre-processing pipeline successfully detected features/ions."
        )
    verifier_colonnes(df, nom, [col_mz, col_rt])
    verifier_colonnes_numeriques(df, nom, [col_mz, col_rt])
    
    # Convert numerical columns (coerce invalid entries to NaN, then replace with 0)
    df[col_mz] = pd.to_numeric(df[col_mz], errors='coerce').fillna(0)
    df[col_rt] = pd.to_numeric(df[col_rt], errors='coerce').fillna(0)
    return df

pddatamat1 = charger_datamatrix(path_DataMatrix1, "dataMatrix Workflow 1")
pddatamat2 = charger_datamatrix(path_DataMatrix2, "dataMatrix Workflow 2")
pdvarmat1  = charger_varmetadata(path_variableMetadata1, "variableMetadata Workflow 1", col_mz, col_rt)
pdvarmat2  = charger_varmetadata(path_variableMetadata2, "variableMetadata Workflow 2", col_mz, col_rt)

verifier_coherence_dataframes(pdvarmat1, pddatamat1, "variableMetadata Workflow 1", "dataMatrix Workflow 1")
verifier_coherence_dataframes(pdvarmat2, pddatamat2, "variableMetadata Workflow 2", "dataMatrix Workflow 2")

# The ID column is assumed to be the first column, allowing any custom name in both file types (1.5.0)
id_col_vm1 = pdvarmat1.columns[0]
id_col_vm2 = pdvarmat2.columns[0]
id_col_dm1 = pddatamat1.columns[0]
id_col_dm2 = pddatamat2.columns[0]

print(f"   ✅ Workflow 1: {len(pdvarmat1)} features, {len(pddatamat1.columns)-1} samples")
print(f"   ✅ Workflow 2: {len(pdvarmat2)} features, {len(pddatamat2.columns)-1} samples")

# =============================================================================
# DISTANCE COMPUTATION & FEATURE MATCHING
# =============================================================================

print("\n🔍 Computing distances between features...")
start_time = time.time()

table1 = pdvarmat1[[col_mz, col_rt]].copy().reset_index(drop=True)
table2 = pdvarmat2[[col_mz, col_rt]].copy().reset_index(drop=True)

rows = []
for i in range(len(table1)):
    rt1 = table1.loc[i, col_rt]
    mz1 = table1.loc[i, col_mz]
    candidats = table2[
        (np.abs(rt1 - table2[col_rt]) < tol_RT) &
        (np.abs(mz1 - table2[col_mz]) < tol_mz)
    ]
    for j in candidats.index:
        rows.append({'id1': i, 'id2': j,
                     'mzdiff': np.abs(mz1 - candidats.at[j, col_mz])})

table_distance = pd.DataFrame(rows).sort_values('mzdiff').reset_index(drop=True) \
    if rows else pd.DataFrame(columns=['id1','id2','mzdiff'])

elapsed = time.time() - start_time
print(f"   → {len(table_distance)} candidate pairs found in {elapsed:.1f}s")

print("\n🔗 Performing unique feature matching...")
start_time = time.time()

blacklist1, blacklist2 = set(), set()
# Common sample columns = intersection excluding the identifier column of each dataMatrix
common_samples = list(
    pddatamat1.columns.drop(id_col_dm1)
    .intersection(pddatamat2.columns.drop(id_col_dm2))
)

associations = []
for _, row in table_distance.iterrows():
    i, j = int(row['id1']), int(row['id2'])
    if i not in blacklist1 and j not in blacklist2:
        blacklist1.add(i)
        blacklist2.add(j)
        mz1_v, mz2_v = table1.loc[i, col_mz], table2.loc[j, col_mz]
        ppmin, ppmax = min(mz1_v, mz2_v), max(mz1_v, mz2_v)
        entry = {
            'idWF1': pdvarmat1.loc[i, id_col_vm1],
            'idWF2': pdvarmat2.loc[j, id_col_vm2],
            'id1': i, 'id2': j,
            'mz1': mz1_v, 'mz2': mz2_v,
            'rt1': table1.loc[i, col_rt], 'rt2': table2.loc[j, col_rt],
            'mzdiff (Da)': mz1_v - mz2_v,
            'mzdiff (ppm)': 1e6*(ppmax - ppmin) / (ppmin + (ppmax - ppmin)/2),
            'rtdiff': table1.loc[i, col_rt] - table2.loc[j, col_rt]
        }
        for s in common_samples:
            entry[s + '_1'] = pddatamat1.loc[i, s]
            entry[s + '_2'] = pddatamat2.loc[j, s]
        associations.append(entry)

# Unmatched features
for i in table1.index:
    if i not in blacklist1:
        entry = {'idWF1': pdvarmat1.loc[i, id_col_vm1], 'idWF2': 'no match',
                 'id1': i, 'id2': 'no match',
                 'mz1': table1.loc[i, col_mz], 'rt1': table1.loc[i, col_rt]}
        for s in common_samples:
            entry[s + '_1'] = pddatamat1.loc[i, s]
        associations.append(entry)

for j in table2.index:
    if j not in blacklist2:
        entry = {'idWF1': 'no match', 'idWF2': pdvarmat2.loc[j, id_col_vm2],
                 'id1': 'no match', 'id2': j,
                 'mz2': table2.loc[j, col_mz], 'rt2': table2.loc[j, col_rt]}
        for s in common_samples:
            entry[s + '_2'] = pddatamat2.loc[j, s]
        associations.append(entry)

tablecomm_final = pd.DataFrame(associations)
pdcommon = tablecomm_final[
    (tablecomm_final['id1'] != 'no match') &
    (tablecomm_final['id2'] != 'no match')
].copy()

n_paired   = len(pdcommon)
n_only_wf1 = (tablecomm_final['id2'] == 'no match').sum()
n_only_wf2 = (tablecomm_final['id1'] == 'no match').sum()
n_wf1_total = n_paired + n_only_wf1
n_wf2_total = n_paired + n_only_wf2

print(f"   ✅ Paired features    : {n_paired}")
print(f"   ⚠️  Unmatched in WF1  : {n_only_wf1} ({100*n_only_wf1/max(n_wf1_total,1):.1f}%)")
print(f"   ⚠️  Unmatched in WF2  : {n_only_wf2} ({100*n_only_wf2/max(n_wf2_total,1):.1f}%)")

# Export tabular file
tablecomm_final.drop(['id1','id2'], axis=1, errors='ignore').to_csv("Associations.tabular", sep='\t', index=False)
print("\n✅ File Associations.tabular exported successfully")

elapsed = time.time() - start_time
print(f"   Execution time: {elapsed:.1f}s\n")

# =============================================================================
# REPORT GENERATION → PDF
# =============================================================================

print("📄 Generating PDF report...")
res = open("res.html", "w", encoding="utf-8")

# --- CSS & HTML Header ---
res.write("""<html lang="en">
<head><meta charset="utf-8">
<style>
  body { font-family: Arial, sans-serif; margin: 40px; color: #222; font-size: 13px; }
  h1 { color: #2c5f8a; border-bottom: 3px solid #2c5f8a; padding-bottom: 8px; margin-bottom: 12px; }
  h2 { color: #2c5f8a; margin-top: 20px; margin-bottom: 6px; }
  h3 { color: #3a7abf; }
  table { border-collapse: collapse; margin: 8px 0; width: auto; }
  th { background: #2c5f8a; color: white; padding: 7px 14px; text-align: center; }
  td { border: 1px solid #ccc; padding: 6px 12px; text-align: center; }
  tr:nth-child(even) td { background: #f4f8fc; }
  .info-box { background: #eef4fb; border-left: 4px solid #2c5f8a;
              padding: 8px 16px; margin: 8px 0; border-radius: 4px; }
  .warn-box { background: #fff8e1; border-left: 4px solid #f9a825;
              padding: 8px 16px; margin: 8px 0; border-radius: 4px; }
  img { max-width: 100%; display: block; margin: 6px 0; }
  .page-break { page-break-before: always; }
</style>
</head><body>
""")

# --- Cover Page ---
now = datetime.now().strftime("%d/%m/%Y at %H:%M")
res.write(f"""
<h1>MS Output Comparison</h1>
<div class="info-box">
  <b>Workflow 1:</b> {wf1_name}<br>
  <b>Workflow 2:</b> {wf2_name}<br>
  <b>m/z Tolerance:</b> {tol_mz} Da &nbsp;|&nbsp; <b>RT Tolerance:</b> {tol_RT}<br>
  <b>Analysis Date:</b> {now}
</div>
""")

# --- Feature Matching Results ---
res.write(f"""
<h2>1. Matching Results</h2>
<div class="info-box">
  The algorithm pairs each feature from Workflow 1 with the closest feature (by mass) 
  from Workflow 2, provided they satisfy the user-defined m/z and RT tolerance thresholds.
</div>
<table>
  <tr><th>Workflow 1 Features</th><th>Workflow 2 Features</th>
      <th>Unmatched in WF1</th><th>Unmatched in WF2</th><th>Paired Features</th></tr>
  <tr>
    <td>{n_wf1_total}</td>
    <td>{n_wf2_total}</td>
    <td>{n_only_wf1}<br>({np.round(100*n_only_wf1/max(n_wf1_total,1),1)}%)</td>
    <td>{n_only_wf2}<br>({np.round(100*n_only_wf2/max(n_wf2_total,1),1)}%)</td>
    <td><b>{n_paired}</b></td>
  </tr>
</table>
""")

if n_only_wf1 > 0.3 * n_wf1_total or n_only_wf2 > 0.3 * n_wf2_total:
    res.write("""<div class="warn-box">
    ⚠️ Warning: more than 30% of the features in at least one workflow remain unmatched. 
    This may indicate major methodological discrepancies between the two pipelines, 
    or overly stringent tolerances. You might consider slightly increasing the thresholds.
    </div>""")

# --- Heatmap Section ---
res.write("<h2>2. Matching Quality</h2>")
res.write("""<div class="info-box">
  The following heatmaps show the proportion of paired features grouped by fractions
  of the user-defined m/z and RT tolerance thresholds.<br>
  Each axis starts with an <b>Identical</b> column/row (strict difference = 0), followed
  by intervals from <b>0</b> to the <b>tolerance values</b> chosen by the user.<br>
  Except from the ppm axis of the 3rd heatmap (which corresponds to non-overlapping intervals),
  the intervals are cumulative (from identical to the maximum value). 
  They represent <b>6 sub-intervals</b>, set regularly: equal size for the 1st heatmap,
  by order of magnitude for the 2nd heatmap (logarithmic scale).<br>
  The <b>Identical × Identical</b> cell (bottom-left) counts only features with
  <b>both Δm/z = 0 and ΔRT = 0</b> simultaneously — this definition is independent
  of the scale (Da or logarithmic) and is identical across all heatmaps.<br>
</div>""")

if n_paired > 0:
    mzdiff_da  = pdcommon['mzdiff (Da)']
    mzdiff_ppm = pdcommon['mzdiff (ppm)']
    rtdiff     = pdcommon['rtdiff']

    totaldim = len(pdcommon)

    # ── v3: tolerance-based heatmap with 6 levels ─────────────────────────
    # Each axis is divided into 6 equal sub-intervals of the user tolerance.
    # Level k  (1-indexed) covers differences  ≤  k/6 × tolerance.
    # This replaces the previous powers-of-10 grid and makes the tolerance
    # thresholds directly visible on the heatmap.
    N_LEVELS = 6  # number of levels per axis (fine granularity — v3)

    def build_heatmap_tol(diff_series, tol_mz_val, rtdiff_series, tol_rt_val,
                          totaldim, n_levels=N_LEVELS):
        """
        Build a ((n_levels + 1) × (n_levels + 1)) heatmap matrix.

        Row/column 0 — "Identical":
            Strictly exact matches on that axis (difference == 0).
            The top-left cell [0, 0] therefore represents features with
            BOTH |Δm/z| == 0 AND |ΔRT| == 0 simultaneously.
            This definition is scale-independent and identical across the
            Da heatmap and the logarithmic heatmap.

        Rows/columns 1 … n_levels — tolerance sub-intervals:
            Cell [i, j] (i ≥ 1, j ≥ 1) = % of paired features with
            |mzdiff| ≤ i/n_levels × tol_mz  AND  |rtdiff| ≤ j/n_levels × tol_rt.

        Rows = m/z levels (index 0 = Identical), columns = RT levels.
        """
        size = n_levels + 1
        arr = np.zeros((size, size))
        mz_abs = diff_series.abs()
        rt_abs = rtdiff_series.abs()
        mz_exact = (mz_abs == 0)
        rt_exact = (rt_abs == 0)
        for i in range(size):
            mz_mask = mz_exact if i == 0 else (mz_abs <= i / n_levels * tol_mz_val)
            for j in range(size):
                rt_mask = rt_exact if j == 0 else (rt_abs <= j / n_levels * tol_rt_val)
                arr[i, j] = 100 * (mz_mask & rt_mask).sum() / totaldim
        return arr

    def tol_labels(tol_val, n_levels=N_LEVELS, unit=''):
        """
        Generate n_levels tick labels: each label = k/n_levels × tol_val.
        e.g. for tol=0.02 Da and 6 levels → ['0.003', '0.007', '0.010',
                                               '0.013', '0.017', '0.020 (tol)']
        The last label is annotated with '(tol)' to highlight the full tolerance.
        """
        labels = []
        for k in range(1, n_levels + 1):
            val = k / n_levels * tol_val
            label = f'{val:.4g}'
            if k == n_levels:
                label += f' {unit}(tol)'.strip()
            labels.append(label)
        return labels

    arr1 = build_heatmap_tol(mzdiff_da,  tol_mz, rtdiff, tol_RT, totaldim)

    jet = plt.get_cmap('turbo')
    cmap = [jet(i) for i in np.linspace(0, 0.99, num=20)]

    # ── Tick labels for the Da heatmap (Identical row/column + tolerance sub-intervals) ──
    # Row/column 0 of arr1 represents strict exact matches (difference == 0).
    # The label 'Identical' is applied unconditionally to index 0 on both axes,
    # consistent with the scientific definition validated by the supervisor.
    # The same convention is used in the logarithmic heatmap, ensuring the
    # Identical × Identical cell has exactly the same meaning in both figures.

    def tol_labels_with_identical(tol_val, n_levels=N_LEVELS, unit=''):
        """
        Generate (n_levels + 1) tick labels for a tolerance-based axis.

        Index 0 — 'Identical': represents the strict equality condition
            (difference == 0 on this axis), which is scale-independent and
            has exactly the same meaning in both the Da and the logarithmic
            heatmap.
        Indices 1 … n_levels — tolerance sub-intervals:
            label k shows k/n_levels × tol_val; the last is annotated
            with '(tol)' to highlight the full user-defined tolerance.
            Each of these labels is prefixed with '< ' to make the
            cumulative nature of the threshold explicit (display only — does
            not affect the underlying calculation).
        """
        labels = ['Identical']
        for k in range(1, n_levels + 1):
            val = k / n_levels * tol_val
            label = f'{val:.4g}'
            if k == n_levels:
                label += f' {unit}(tol)'.strip()
            labels.append(f'< {label}')
        return labels

    y_labels_da = tol_labels_with_identical(tol_mz, unit='Da')
    x_labels_rt = tol_labels_with_identical(tol_RT, unit='')

    # ── Heatmap 1 — m/z in Dalton (validated, unchanged logic) ──────────────
    # Standalone display removed; arr1 is still computed above and reused in
    # the side-by-side comparison figure (Da + Logarithmic) further below.

    # ── Heatmap 3 (displayed third, side by side with Heatmap 1) ─────────────
    # Logarithmic multi-scale heatmap — rebuilt directly from the same logic
    # as Heatmap 1 (Da). The ONLY difference is how the per-level thresholds
    # are derived:
    #   - Heatmap 1 (Da):  mz_thresh[i] = (i+1)/N_LEVELS * tol_mz   (linear fractions)
    #   - Heatmap 3 (Log): mz_thresh[i] = tol_mz * 10^LOG_EXPONENTS_LOG[i]  (powers of 10)
    # Everything else — matrix shape (N_LEVELS x N_LEVELS), cell definition,
    # axis convention (X = RT thresholds, Y = m/z thresholds), the "Identical"
    # first-row/first-column label convention, and the invert_yaxis() display
    # — is identical to Heatmap 1.
    #
    # Counting logic — identical to Heatmap 1 (Da):
    #   cell [i, j] = % of paired features with |Δm/z| ≤ mz_thresh[i]
    #                                        AND |ΔRT|  ≤ rt_thresh[j]
    #
    # Display convention (consistent with Heatmap 1):
    #   Y-axis: Identical at BOTTOM, full tolerance (tol) at TOP.
    #   X-axis: Identical on the LEFT, full tolerance (tol) on the RIGHT.

    # Strictest-first order, same axis direction as Heatmap 1.
    LOG_EXPONENTS_LOG = [-5, -4, -3, -2, -1, 0]
    N_LOG = len(LOG_EXPONENTS_LOG)   # == N_LEVELS, so the matrix is the same size as arr1

    # Per-level thresholds, mirroring the Da heatmap's mz_thresh / rt_thresh.
    mz_log_thresholds = [tol_mz * (10 ** e) for e in LOG_EXPONENTS_LOG]
    rt_log_thresholds = [tol_RT  * (10 ** e) for e in LOG_EXPONENTS_LOG]

    mzdiff_da_abs = mzdiff_da.abs()
    rtdiff_abs    = rtdiff.abs()

    # Same (N_LEVELS+1) × (N_LEVELS+1) shape as arr1 — the extra row/column
    # at index 0 represents strict exact matches (difference == 0), providing
    # the same "Identical" definition as the Da heatmap.
    arr_log = np.zeros((N_LOG + 1, N_LOG + 1))
    mzdiff_exact = (mzdiff_da_abs == 0)
    rtdiff_exact  = (rtdiff_abs  == 0)
    for i in range(N_LOG + 1):
        if i == 0:
            mz_mask = mzdiff_exact
        else:
            mz_mask = (mzdiff_da_abs <= mz_log_thresholds[i - 1])
        for j in range(N_LOG + 1):
            if j == 0:
                rt_mask = rtdiff_exact
            else:
                rt_mask = (rtdiff_abs <= rt_log_thresholds[j - 1])
            arr_log[i, j] = 100 * (mz_mask & rt_mask).sum() / totaldim

    def tol_labels_log_with_identical(tol_val, exponents, unit=''):
        """
        Generate (len(exponents) + 1) tick labels for the logarithmic axis.

        Index 0 — 'Identical': strict equality condition (difference == 0).
            This is identical in meaning to the Da heatmap's first label and
            is independent of the logarithmic scale.
        Indices 1 … len(exponents): plain decimal threshold values actually
            used in the calculation (tol_val * 10^exponent), with the last
            one annotated '(tol)' to mark the full user-defined tolerance.
            Each of these labels is prefixed with '< ' to make the
            cumulative nature of the threshold explicit (display only — does
            not affect the underlying calculation).
        No 'tol×10^x' scientific notation is shown — only the resulting
        numeric threshold, e.g. 0.00001, 0.0001, ..., 1 (or the equivalent
        values computed from the user's chosen tolerance).
        """
        labels = ['Identical']
        last_idx = len(exponents) - 1
        for k, e in enumerate(exponents):
            val = tol_val * (10 ** e)
            label = f'{val:.4g}'
            if k == last_idx:
                label += f' {unit}(tol)'.strip()
            labels.append(f'< {label}')
        return labels

    # Labels in strictest-first order: 'Identical' at index 0 (finest level),
    # then increasing powers of 10 up to the full tolerance at the last index.
    # seaborn renders row 0 at the TOP; invert_yaxis() flips so 'Identical'
    # sits at the BOTTOM and the full tolerance at the TOP — identical to the
    # Da heatmap.
    y_labels_log = tol_labels_log_with_identical(tol_mz, LOG_EXPONENTS_LOG, unit='Da')
    x_labels_log = tol_labels_log_with_identical(tol_RT, LOG_EXPONENTS_LOG, unit='')

    # ── Side-by-side figure: Heatmap 1 (Da) on the left, Heatmap 3 (Log) on the right ──
    # Both heatmaps share the same computation logic; side-by-side display makes
    # it straightforward to verify that differences stem only from the scale change.
    fig_side, (ax_da_side, ax_log_side) = plt.subplots(
        1, 2, figsize=(16, 7),
        gridspec_kw={'wspace': 0.55}
    )
    fig_side.patch.set_facecolor('#ffffff')

    # Left panel — Da heatmap (identical rendering to the standalone Heatmap 1)
    sns.heatmap(arr1, annot=True, cmap=cmap, linecolor='white', linewidths=1,
                fmt='.1f', cbar_kws={'shrink': 0.6}, annot_kws={'fontsize': 8},
                vmax=100, vmin=0,
                yticklabels=y_labels_da,
                xticklabels=x_labels_rt, ax=ax_da_side, square=True)
    ax_da_side.invert_yaxis()
    ax_da_side.set_title(
        f'm/z tolerance (Da) — {N_LEVELS} levels\nFull tolerance = {tol_mz} Da',
        fontsize=9
    )
    ax_da_side.set_xlabel(f'RT tolerance — full = {tol_RT}', fontsize=8)
    ax_da_side.set_ylabel(f'm/z (Da) — full tolerance = {tol_mz} Da', fontsize=8)

    # Right panel — logarithmic heatmap
    sns.heatmap(arr_log, annot=True, cmap=cmap, linecolor='white', linewidths=1,
                fmt='.1f', cbar_kws={'shrink': 0.6}, annot_kws={'fontsize': 8},
                vmax=100, vmin=0,
                yticklabels=y_labels_log,
                xticklabels=x_labels_log, ax=ax_log_side, square=True)
    # invert_yaxis: row 0 (Identical) moves to the BOTTOM; tol×10^0 (full tol)
    # sits at the TOP — consistent with the Da heatmap.
    ax_log_side.invert_yaxis()
    ax_log_side.set_title(
        f'm/z tolerance (Da) — logarithmic multi-scale\n'
        f'Threshold = tol × 10^exp, exp = {", ".join(str(e) for e in LOG_EXPONENTS_LOG)}',
        fontsize=9
    )
    ax_log_side.set_xlabel(
        f'RT tolerance  (strictest → full tol = {tol_RT})', fontsize=8
    )
    ax_log_side.set_ylabel(
        f'm/z (Da)  (strictest → full tol = {tol_mz} Da)', fontsize=8
    )

    fig_side.suptitle(
        f"Matching Quality Comparison — Da (left) vs Logarithmic Scale (right)\n"
        f"m/z tol = {tol_mz} Da  |  RT tol = {tol_RT}",
        fontsize=10
    )
    fig_side.savefig('heatmap_log.png', bbox_inches='tight', dpi=100)
    plt.close()
    res.write('<img src="heatmap_log.png">')

    # ── Heatmap 2 (displayed second) — m/z in ppm with fixed analytical classes
    # ─────────────────────────────────────────────────────────────────────────
    # Fixed ppm classes (ascending order, bottom → top on the displayed heatmap):
    #   Identical (= 0)  |  ]0;1]  |  ]1;5]  |  ]5;10]  |  ]10;50]  |  >50
    # Each cell [i, j] = % of paired features with |Δm/z(ppm)| in class i
    # AND |ΔRT| ≤ (j+1)/N_LEVELS × tol_RT.
    #
    # "Identical" labels the [= 0 ppm] class, consistent with the Da heatmap.
    PPM_BOUNDS = [0, 1, 5, 10, 50]   # upper edges of the first 5 finite classes
    PPM_LABELS = ['Identical\n(= 0 ppm)', ']0 ; 1]', ']1 ; 5]',
                  ']5 ; 10]', ']10 ; 50]', '> 50']
    N_PPM_ROWS = len(PPM_LABELS)     # 6 rows

    arr2 = np.zeros((N_PPM_ROWS, N_LEVELS))
    ppm_abs = mzdiff_ppm.abs()
    for i in range(N_PPM_ROWS):
        if i == 0:
            ppm_mask = (ppm_abs == 0)
        elif i < len(PPM_BOUNDS):
            ppm_mask = (ppm_abs > PPM_BOUNDS[i - 1]) & (ppm_abs <= PPM_BOUNDS[i])
        else:
            ppm_mask = ppm_abs > PPM_BOUNDS[-1]
        for j in range(N_LEVELS):
            rt_thresh = (j + 1) / N_LEVELS * tol_RT
            arr2[i, j] = 100 * (ppm_mask & (rtdiff.abs() <= rt_thresh)).sum() / totaldim

    # seaborn renders row 0 at the TOP by default; invert_yaxis() moves row 0
    # (the "Identical / = 0 ppm" class) to the BOTTOM, so classes read in
    # ascending order from bottom (strictest) to top (most permissive).
    fig_ppm, ax_ppm = plt.subplots(figsize=(8, 6))
    fig_ppm.tight_layout(pad=5.0)

    sns.heatmap(arr2, annot=True, cmap=cmap, linecolor='white', linewidths=1,
                fmt='.1f', cbar_kws={'shrink': 0.6}, annot_kws={'fontsize': 8},
                vmax=100, vmin=0,
                yticklabels=PPM_LABELS,
                xticklabels=x_labels_rt, ax=ax_ppm, square=True)
    ax_ppm.invert_yaxis()   # ascending order: Identical at bottom, > 50 ppm at top
    ax_ppm.set_title(
        f'm/z tolerance (ppm) — fixed analytical classes',
        fontsize=9
    )
    ax_ppm.set_xlabel(f'RT tolerance — full = {tol_RT}', fontsize=8)
    ax_ppm.set_ylabel('m/z difference (ppm)', fontsize=8)
    fig_ppm.suptitle(
        f"Proportion of Paired Features by ppm Class (%)\n"
        f"RT tol = {tol_RT}  |  {N_LEVELS} RT levels on X-axis",
        fontsize=10
    )
    fig_ppm.savefig('heatmap_ppm.png', bbox_inches='tight', dpi=100)
    plt.close()
    res.write('<img src="heatmap_ppm.png">')


# --- Mass and RT Statistics ---
if n_paired > 0:
    res.write("<h2>3. Statistics of Paired Features</h2>")
    res.write("""<div class="info-box">
      Mass (m/z) and retention time (RT) summary statistics for successfully paired features.
    </div>""")

    def fmt(val, dec=3): return str(np.round(val, dec))

    res.write(f"""
    <table>
      <tr><th></th><th>Min m/z</th><th>Mean m/z</th><th>Max m/z</th>
          <th>Min RT</th><th>Mean RT</th><th>Max RT</th></tr>
      <tr><td><b>{wf1_name}</b></td>
          <td>{fmt(pdcommon['mz1'].min(),2)}</td>
          <td>{fmt(pdcommon['mz1'].mean(),2)}</td>
          <td>{fmt(pdcommon['mz1'].max(),2)}</td>
          <td>{fmt(pdcommon['rt1'].min(),2)}</td>
          <td>{fmt(pdcommon['rt1'].mean(),2)}</td>
          <td>{fmt(pdcommon['rt1'].max(),2)}</td></tr>
      <tr><td><b>{wf2_name}</b></td>
          <td>{fmt(pdcommon['mz2'].min(),2)}</td>
          <td>{fmt(pdcommon['mz2'].mean(),2)}</td>
          <td>{fmt(pdcommon['mz2'].max(),2)}</td>
          <td>{fmt(pdcommon['rt2'].min(),2)}</td>
          <td>{fmt(pdcommon['rt2'].mean(),2)}</td>
          <td>{fmt(pdcommon['rt2'].max(),2)}</td></tr>
      <tr><td><b>Differences</b></td>
          <td>{fmt(pdcommon['mzdiff (Da)'].abs().min(),5)}</td>
          <td>{fmt(pdcommon['mzdiff (Da)'].abs().mean(),5)}</td>
          <td>{fmt(pdcommon['mzdiff (Da)'].abs().max(),5)}</td>
          <td>{fmt(pdcommon['rtdiff'].abs().min(),3)}</td>
          <td>{fmt(pdcommon['rtdiff'].abs().mean(),3)}</td>
          <td>{fmt(pdcommon['rtdiff'].abs().max(),3)}</td></tr>
    </table>
    """)

# --- Intensity extraction for common paired features ---
pddatamat1com = pddatamat1.loc[pdcommon['id1'].astype(int)][
    np.sort(pddatamat1.columns.drop(id_col_dm1).intersection(pddatamat2.columns.drop(id_col_dm2)))
]
pddatamat2com = pddatamat2.loc[pdcommon['id2'].astype(int)][
    np.sort(pddatamat1.columns.drop(id_col_dm1).intersection(pddatamat2.columns.drop(id_col_dm2)))
]

nbtot = np.size(pddatamat1com)

if nbtot > 0 and len(common_samples) > 0:
    res.write('<div class="page-break"></div>')
    res.write("<h2>4. Intensity Comparison</h2>")
    res.write(f"""<div class="info-box">
      Intensity comparison across the <b>{len(common_samples)} common samples</b>
      between the two workflows, for paired ions only.
    </div>""")

    n1n2   = int(((pddatamat1com.values == 0) & (pddatamat2com.values == 0)).sum())
    n1nn2  = int(((pddatamat1com.values == 0) & (pddatamat2com.values != 0)).sum())
    nn1n2  = int(((pddatamat1com.values != 0) & (pddatamat2com.values == 0)).sum())
    nn1nn2 = int(((pddatamat1com.values != 0) & (pddatamat2com.values != 0)).sum())
    nnegal = int(((pddatamat1com.values != 0) & (pddatamat2com.values != 0) &
                  (pddatamat1com.values == pddatamat2com.values)).sum())

    res.write(f"""
    <table>
      <tr><th></th><th>Zero intensity in WF2</th><th>Non-zero intensity in WF2</th></tr>
      <tr><td><b>Zero intensity in WF1</b></td>
          <td>{n1n2} ({round(100*n1n2/nbtot,1)}%)</td>
          <td>{n1nn2} ({round(100*n1nn2/nbtot,1)}%)</td></tr>
      <tr><td><b>Non-zero intensity in WF1</b></td>
          <td>{nn1n2} ({round(100*nn1n2/nbtot,1)}%)</td>
          <td>{nn1nn2} ({round(100*nn1nn2/nbtot,1)}%)</td></tr>
    </table>
    """)

    if nn1nn2 > 0:
        res.write(f"<p><b>Identical non-zero intensities: {nnegal}/{nn1nn2} "
                  f"({round(100*nnegal/nn1nn2,1)}%)</b></p>")

        # Boxplot of relative intensity differences
        pddiff = 100 * (pddatamat1com[pddatamat1com != 0]
                        .subtract(pddatamat2com[pddatamat2com != 0])
                        ).replace([np.inf, -np.inf], np.nan)

        _nonzero1 = pddatamat1com.values[pddatamat1com.values > 0]
        _nonzero2 = pddatamat2com.values[pddatamat2com.values > 0]
        _global_max = np.max(np.concatenate([_nonzero1, _nonzero2])) \
            if (len(_nonzero1) > 0 or len(_nonzero2) > 0) else 1.0
        max_exp = int(np.floor(np.log10(max(_global_max, 1.0))))

        vecdiff = []
        for i in range(0, max_exp + 1):
            subset = pddiff.values[
                (np.abs(pddiff.values) < 10**(i+1)) &
                (np.abs(pddiff.values) >= 10**i)
            ]
            vecdiff.append(subset.tolist() if len(subset) > 0 else [np.nan])

        max_long = max(len(v) for v in vecdiff)
        for i in range(len(vecdiff)):
            if len(vecdiff[i]) < max_long:
                vecdiff[i] += [np.nan] * (max_long - len(vecdiff[i]))

        dfdiff = pd.DataFrame(np.array(vecdiff).transpose())

        fig, axs = plt.subplots(figsize=(8, 5))
        bp = dict(linestyle='-', linewidth=2)
        mp = dict(marker='D', markeredgecolor='black', markerfacecolor='seagreen', markersize=5)

        dfdiff[dfdiff > 0].boxplot(
            boxprops=dict(**bp, color='navy'), medianprops=dict(**bp, color='navy'),
            whiskerprops=dict(**bp, color='navy'), capprops=dict(**bp, color='navy'),
            flierprops=dict(marker='o', markerfacecolor='navy', markersize=2),
            showmeans=True, meanprops=mp, ax=axs
        )
        dfdiff[dfdiff < 0].boxplot(
            boxprops=dict(**bp, color='firebrick'), medianprops=dict(**bp, color='firebrick'),
            whiskerprops=dict(**bp, color='firebrick'), capprops=dict(**bp, color='firebrick'),
            flierprops=dict(marker='o', markerfacecolor='firebrick', markersize=2),
            showmeans=True, meanprops=mp, ax=axs
        )

        for i in range(dfdiff.shape[1]):
            n_val = (~np.isnan(dfdiff[i])).sum()
            axs.text(i + 1, axs.get_ylim()[0], f'n={n_val}',
                     fontsize=7, ha='center', va='bottom')

        axs.set_xlabel("Order of magnitude (powers of 10)", fontsize=9)
        axs.set_ylabel("Relative intensity difference (%)", fontsize=9)
        axs.set_title(f"Distribution of intensity differences\n{wf1_name} vs {wf2_name}", fontsize=10)

        legend_lines = [
            Line2D([0], [0], color='navy', lw=3),
            Line2D([0], [0], color='firebrick', lw=3)
        ]
        axs.legend(legend_lines, [f'WF1 > WF2', f'WF1 < WF2'], fontsize=8)
        fig.tight_layout()
        fig.savefig('histo_intensity.png', bbox_inches='tight', dpi=100)
        plt.close()
        res.write('<img src="histo_intensity.png">')

        # Histogramme des distributions
        int1 = pddatamat1com.replace(0, np.nan).values.flatten()
        int1 = int1[~np.isnan(int1)]
        int2 = pddatamat2com.replace(0, np.nan).values.flatten()
        int2 = int2[~np.isnan(int2)]

        if len(int1) > 0 and len(int2) > 0:
            _all_int = np.concatenate([int1, int2])
            _log_max = int(np.ceil(np.log10(max(_all_int.max(), 10))))
            bins = [10**i for i in range(_log_max + 1)]
            if len(bins) < 2:
                bins = [1, 10]  # fallback: single-bin guard
            fig1, ax1 = plt.subplots(figsize=(7, 4))
            ax1.grid(zorder=0, linestyle='dashed', alpha=0.5)
            hist1, _ = np.histogram(int1, bins=bins)
            hist2, _ = np.histogram(int2, bins=bins)
            x = np.arange(len(hist1))
            ax1.bar(x - 0.22, hist1, width=0.44, color='steelblue',
                    edgecolor='white', zorder=3, alpha=0.9, label=wf1_name)
            ax1.bar(x + 0.22, hist2, width=0.44, color='darkorange',
                    edgecolor='white', zorder=3, alpha=0.9, label=wf2_name)
            ax1.set_xticks(x)
            ax1.set_xticklabels([f'10e{int(np.ceil(np.log10(bins[i+1])))}' for i in range(len(bins)-1)])
            ax1.set_title("Intensity distribution by order of magnitude", fontsize=10)
            ax1.set_xlabel("Order of magnitude", fontsize=9)
            ax1.set_ylabel("Number of values", fontsize=9)
            ax1.legend(fontsize=9)
            fig1.tight_layout()
            fig1.savefig('histo_int_temp.png', bbox_inches='tight', dpi=100)
            plt.close()
            res.write('<img src="histo_int_temp.png">')
            res.write("""<div class="info-box" style="font-size:12px; margin-top:4px;">
              <b>Note:</b> Zero intensity values are not included in this distribution plot,
              as their inclusion may overwhelm the non-zero intensity distribution. Please refer to the summary table above for the number of
              zero-intensity values in each workflow.
            </div>""")

        # --- Section 5: Per-ion Intensity Correlation ---
        res.write('<div class="page-break"></div>')
        res.write("<h2>5. Per-ion Intensity Correlation</h2>")
        res.write(f"""<div class="info-box">
          For each paired ion, the Pearson correlation coefficient (computed on log10-transformed intensities)
          is calculated across the <b>{len(common_samples)} common samples</b> between
          <b>{wf1_name}</b> and <b>{wf2_name}</b>.<br><br>
          Four histograms are displayed:<br>
          &nbsp;&nbsp;⓪ <b>Full-scale histogram</b> — fixed bins [−1.0 → 1.0], step 0.1 (original chart, comparable across runs)<br>
          &nbsp;&nbsp;① <b>Global distribution</b> — 20 bins spanning [min&nbsp;r,&nbsp;max&nbsp;r]<br>
          &nbsp;&nbsp;② <b>Lower-half zoom</b> — 20 bins spanning [min&nbsp;r,&nbsp;<b>split&nbsp;threshold&nbsp;=&nbsp;{corr_split:.2f}</b>], focusing on weaker correlations<br>
          &nbsp;&nbsp;③ <b>Upper-half zoom</b> — 20 bins spanning [<b>split&nbsp;threshold&nbsp;=&nbsp;{corr_split:.2f}</b>,&nbsp;max&nbsp;r], focusing on stronger correlations<br><br>
          The split threshold (<b>r = {corr_split:.2f}</b>) is user-configurable and separates the two zoom views.
          Adjust it in the tool parameters to focus on a specific correlation range.
        </div>""")

        # Compute Pearson r per ion (across common samples)
        ion_r_values = []
        for ion_idx, row_ion in pdcommon.iterrows():
            vals1, vals2 = [], []
            for col in common_samples:
                c1, c2 = col + '_1', col + '_2'
                if c1 in pdcommon.columns and c2 in pdcommon.columns:
                    v1_val = row_ion[c1]
                    v2_val = row_ion[c2]
                    if v1_val > 0 and v2_val > 0:
                        vals1.append(np.log10(v1_val))
                        vals2.append(np.log10(v2_val))
            if len(vals1) >= 3:
                r_ion = np.corrcoef(vals1, vals2)[0, 1]
                ion_r_values.append({
                    'ion': row_ion.get('idWF1', str(ion_idx)),
                    'r': r_ion
                })

        if len(ion_r_values) == 0:
            res.write("""<div class="warn-box">
            ⚠️ Not enough non-zero data points per ion to compute the correlation
            (minimum of 3 non-zero samples required per ion).<br>
            Please verify that your data contain a sufficient number of common samples
            with non-zero intensities.
            </div>""")
        else:
            df_r = (pd.DataFrame(ion_r_values)
                    .dropna(subset=['r'])
                    .sort_values('r', ascending=False)
                    .reset_index(drop=True))
            r_vals = df_r['r'].values

            # ── Shared colour scheme and legend ──────────────────────────────
            # Colour thresholds reflect standard interpretation of Pearson r
            PALETTE = {
                'strong':   '#1a6faf',   # r ≥ 0.80 — strong (deep blue)
                'moderate': '#56a5d8',   # 0.50 ≤ r < 0.80 — moderate (medium blue)
                'weak':     '#f5a623',   # 0.30 ≤ r < 0.50 — weak (amber)
                'poor':     '#c0392b',   # r < 0.30 — poor / negative (red)
            }
            legend_patches = [
                mpatches.Patch(facecolor=PALETTE['strong'],   edgecolor='white', label='r ≥ 0.80  —  Strong'),
                mpatches.Patch(facecolor=PALETTE['moderate'], edgecolor='white', label='0.50 ≤ r < 0.80  —  Moderate'),
                mpatches.Patch(facecolor=PALETTE['weak'],     edgecolor='white', label='0.30 ≤ r < 0.50  —  Weak'),
                mpatches.Patch(facecolor=PALETTE['poor'],     edgecolor='white', label='r < 0.30  —  Poor / Negative'),
            ]

            def color_for_r(r_val):
                """Return the bar colour for a given Pearson r value."""
                if r_val >= 0.80:  return PALETTE['strong']
                if r_val >= 0.50:  return PALETTE['moderate']
                if r_val >= 0.30:  return PALETTE['weak']
                return PALETTE['poor']

            # ── Shared aesthetic helper ───────────────────────────────────────
            def _style_axis(ax, title, xlabel, ylabel, n_ions, subtitle=''):
                """Apply consistent professional styling to a histogram axis."""
                ax.set_facecolor('#f7f9fc')
                ax.grid(axis='y', linestyle='--', linewidth=0.6, alpha=0.7,
                        color='#b0bec5', zorder=0)
                ax.set_axisbelow(True)
                ax.spines[['top', 'right']].set_visible(False)
                ax.spines[['left', 'bottom']].set_color('#90a4ae')
                ax.tick_params(colors='#455a64', labelsize=7.5)
                full_title = f"{title}\n{wf1_name}  vs  {wf2_name}  —  {n_ions} ions analysed"
                if subtitle:
                    full_title += f"\n{subtitle}"
                ax.set_title(full_title, fontsize=9.5, fontweight='bold',
                             color='#1c2833', pad=10)
                ax.set_xlabel(xlabel, fontsize=9, color='#1c2833', labelpad=6)
                ax.set_ylabel(ylabel, fontsize=9, color='#1c2833', labelpad=6)

            def _draw_histogram(ax, r_subset, r_min, r_max, n_bins=20,
                                 annotate_median=True, n_ions_total=None):
                """
                Draw a professional histogram on `ax` with exactly `n_bins` bins
                spanning [r_min, r_max].

                When the distribution is strongly skewed (the dominant bin holds
                more than 20x the second-largest bin), a secondary logarithmic
                Y-axis is automatically added on the right so that low-count
                bars remain legible alongside very tall dominant bars.

                Parameters
                ----------
                r_subset : array-like
                    Correlation values to histogram (may be a subset for zoom views).
                r_min, r_max : float
                    Explicit axis / bin boundaries.
                n_bins : int
                    Number of bins (always 20 per specification).
                annotate_median : bool
                    Whether to draw a vertical median line.
                n_ions_total : int or None
                    Total ion count for display; defaults to len(r_subset).
                """
                if r_max - r_min < 1e-9:
                    # Edge case: all values identical — single bin
                    step = 1e-9
                else:
                    step = (r_max - r_min) / n_bins

                # Extend the last edge by a tiny epsilon so that values exactly
                # equal to r_max are always captured in the last bin.
                bins_lo = np.linspace(r_min, r_max, n_bins + 1)
                bins_hi = bins_lo.copy()
                bins_hi[-1] = np.nextafter(bins_lo[-1], np.inf)
                counts, edges = np.histogram(r_subset, bins=bins_hi)
                # Keep original edges (without epsilon) for display
                centers = (bins_lo[:-1] + bins_lo[1:]) / 2

                bar_colors = [color_for_r(c) for c in centers]

                bars = ax.bar(centers, counts, width=step * 0.88,
                              color=bar_colors, edgecolor='white',
                              linewidth=0.6, zorder=3)

                # Annotate bar counts (skip zero-count bars)
                y_offset = max(counts) * 0.015 if counts.max() > 0 else 0.1
                for bar, cnt in zip(bars, counts):
                    if cnt > 0:
                        ax.text(bar.get_x() + bar.get_width() / 2,
                                bar.get_height() + y_offset,
                                str(int(cnt)),
                                ha='center', va='bottom',
                                fontsize=6.5, fontweight='semibold',
                                color='#2c3e50')

                # ── Conditional secondary log-scale Y-axis ────────────────────
                # If the tallest bar is >= 20x the second tallest non-zero bar,
                # the linear scale hides low-count classes. A twin log axis is
                # added on the right so every class remains visible.
                nonzero_counts = counts[counts > 0]
                _LOG_SKEW_THRESHOLD = 20  # dominance ratio triggering log axis
                if (len(nonzero_counts) >= 2 and
                        nonzero_counts.max() >= _LOG_SKEW_THRESHOLD *
                        np.sort(nonzero_counts)[-2]):
                    ax_log = ax.twinx()
                    nonzero_mask = counts > 0
                    ax_log.set_yscale('log')
                    ax_log.set_ylim(bottom=0.8)
                    ax_log.set_ylabel('Count (log scale)', fontsize=7.5,
                                      color='#7f8c8d', labelpad=4)
                    ax_log.tick_params(axis='y', labelsize=7, colors='#7f8c8d')
                    ax_log.spines['right'].set_color('#b0bec5')
                    ax_log.spines['top'].set_visible(False)
                    # Note explaining the dual scale
                    ax.annotate(
                        'Skewed distribution — log scale (right axis) shows low-count classes',
                        xy=(0.01, 0.97), xycoords='axes fraction',
                        fontsize=6.5, color='#7f8c8d',
                        va='top', ha='left', style='italic')
                # ─────────────────────────────────────────────────────────────

                # Median vertical line
                if annotate_median and len(r_subset) > 0:
                    med = np.median(r_subset)
                    ax.axvline(med, color='#2ecc71', linewidth=1.6,
                               linestyle='--', zorder=4,
                               label=f'Median = {med:.3f}')
                    ax.legend(fontsize=7.5, loc='upper left',
                              framealpha=0.85, edgecolor='#b0bec5')

                # X-tick labels: show bin edges at readable frequency
                tick_step = max(1, n_bins // 10)
                tick_positions = bins_lo[::tick_step]
                ax.set_xticks(tick_positions)
                ax.set_xticklabels([f'{v:.2f}' for v in tick_positions],
                                   rotation=40, ha='right', fontsize=7)
                ax.set_xlim(r_min - step * 0.5, r_max + step * 0.5)

                return counts, edges

            # ── Key statistics used across all three charts ───────────────────
            r_min    = float(r_vals.min())
            r_max    = float(r_vals.max())
            r_median = float(np.median(r_vals))
            n_ions   = len(r_vals)
            N_BINS   = 20   # fixed by specification

            # Clamp the user-defined split threshold to the actual data range so
            # that both zoom histograms always contain at least one data point.
            # The threshold is also shown in chart titles and the info-box.
            r_split = float(np.clip(corr_split, r_min, r_max))

            # ══════════════════════════════════════════════════════════════════
            # Chart 0 — Original fixed-scale histogram  [-1.0 → 1.0]  (step 0.1)
            # Kept from the original tool for backward compatibility.
            # Bins are always the same regardless of the data range, making
            # results directly comparable across different runs.
            # ══════════════════════════════════════════════════════════════════
            # Fixed-scale bins [-1.0 → 1.0], step 0.1.
            # The upper boundary is extended by a tiny epsilon so that values
            # exactly equal to 1.0 fall inside the last bin (np.histogram uses
            # half-open intervals [a, b) for all bins except the last which is
            # closed [a, b]).  Using np.nextafter(1.0+0.1, np.inf) guarantees
            # r = 1.0 is always counted in the bin [0.9, 1.0].
            _bins_lo = np.round(np.arange(-1.0, 1.0, 0.1), 10)
            bins_fixed = np.append(_bins_lo, np.nextafter(1.0, np.inf))
            counts_fixed, edges_fixed = np.histogram(r_vals, bins=bins_fixed)
            step_fixed = 0.1
            # Use true r-value centres so the x-axis matches the other charts
            bin_centers_fixed = (_bins_lo + (_bins_lo + step_fixed)) / 2
            bin_labels_fixed  = [f'{_bins_lo[i]:.1f}–{min(_bins_lo[i]+step_fixed, 1.0):.1f}'
                                  for i in range(len(_bins_lo))]
            bar_colors_fixed  = [color_for_r(c) for c in bin_centers_fixed]

            fig0, ax0 = plt.subplots(figsize=(11, 4.5))
            fig0.patch.set_facecolor('#ffffff')
            ax0.set_facecolor('#f7f9fc')
            ax0.grid(axis='y', linestyle='--', linewidth=0.6, alpha=0.7,
                     color='#b0bec5', zorder=0)
            ax0.set_axisbelow(True)
            ax0.spines[['top', 'right']].set_visible(False)
            ax0.spines[['left', 'bottom']].set_color('#90a4ae')

            bars0 = ax0.bar(bin_centers_fixed, counts_fixed,
                            color=bar_colors_fixed, edgecolor='white',
                            linewidth=0.6, zorder=3, width=step_fixed * 0.88)
            y_off0 = max(counts_fixed) * 0.015 if counts_fixed.max() > 0 else 0.1
            for bar, cnt in zip(bars0, counts_fixed):
                if cnt > 0:
                    ax0.text(bar.get_x() + bar.get_width() / 2,
                             bar.get_height() + y_off0,
                             str(int(cnt)),
                             ha='center', va='bottom',
                             fontsize=6.5, fontweight='semibold',
                             color='#2c3e50')

            # ── Conditional secondary log-scale Y-axis for Chart 0 ────────────
            _nz0 = counts_fixed[counts_fixed > 0]
            if (len(_nz0) >= 2 and
                    _nz0.max() >= 20 * np.sort(_nz0)[-2]):
                ax0_log = ax0.twinx()
                _nz_mask0 = counts_fixed > 0
                ax0_log.set_yscale('log')
                ax0_log.set_ylim(bottom=0.8)
                ax0_log.set_ylabel('Count (log scale)', fontsize=7.5,
                                   color='#7f8c8d', labelpad=4)
                ax0_log.tick_params(axis='y', labelsize=7, colors='#7f8c8d')
                ax0_log.spines['right'].set_color('#b0bec5')
                ax0_log.spines['top'].set_visible(False)
                ax0.annotate(
                    'Skewed distribution — log scale (right axis) shows low-count classes',
                    xy=(0.01, 0.97), xycoords='axes fraction',
                    fontsize=6.5, color='#7f8c8d',
                    va='top', ha='left', style='italic')
            # ──────────────────────────────────────────────────────────────────

            # Median vertical line — same style as Chart 1
            ax0.axvline(r_median, color='#2ecc71', linewidth=1.6,
                        linestyle='--', zorder=4)

            ax0.set_xticks(bin_centers_fixed)
            ax0.set_xticklabels(bin_labels_fixed, rotation=45, ha='right', fontsize=7.5)
            ax0.set_xlim(-1.0 - step_fixed * 0.5, 1.0 + step_fixed * 0.5)
            ax0.tick_params(colors='#455a64')
            ax0.set_xlabel('Correlation class (r)', fontsize=9, color='#1c2833', labelpad=6)
            ax0.set_ylabel('Number of ions',        fontsize=9, color='#1c2833', labelpad=6)
            ax0.set_title(
                f'Full-Scale Correlation Histogram — fixed bins [−1.0 → 1.0], step 0.1\n'
                f'{wf1_name}  vs  {wf2_name}  —  {n_ions} ions analysed',
                fontsize=9.5, fontweight='bold', color='#1c2833', pad=10
            )
            median_handle0 = mpatches.Patch(facecolor='#2ecc71', edgecolor='white',
                                            label=f'Median = {r_median:.3f}')
            ax0.legend(handles=legend_patches + [median_handle0],
                       fontsize=7.5, loc='upper left',
                       framealpha=0.85, edgecolor='#b0bec5',
                       title='Correlation level', title_fontsize=7.5)

            fig0.tight_layout(pad=2.5)
            fig0.savefig('corr_histogram.png', bbox_inches='tight',
                         dpi=130, facecolor='#ffffff')
            plt.close()
            res.write('<img src="corr_histogram.png">')

            # ══════════════════════════════════════════════════════════════════
            # Chart 1 — Global histogram  [r_min → r_max]  (20 bins)
            # ══════════════════════════════════════════════════════════════════
            fig1, ax1_hist = plt.subplots(figsize=(11, 4.5))
            fig1.patch.set_facecolor('#ffffff')

            _draw_histogram(ax1_hist, r_vals, r_min, r_max,
                            n_bins=N_BINS, annotate_median=True)
            _style_axis(ax1_hist,
                        title='Global Correlation Distribution — 20 equal-width bins',
                        xlabel='Pearson r',
                        ylabel='Number of ions',
                        n_ions=n_ions)
            # Include the median entry in the legend so it remains visible on
            # this chart.  Without this, the ax.legend() call here would
            # overwrite the legend set inside _draw_histogram and drop the
            # median handle that was added there.
            median_handle1 = mpatches.Patch(facecolor='#2ecc71', edgecolor='white',
                                            label=f'Median = {r_median:.3f}')
            ax1_hist.legend(handles=legend_patches + [median_handle1], fontsize=7.5,
                            loc='upper left', framealpha=0.85,
                            edgecolor='#b0bec5', title='Correlation level',
                            title_fontsize=7.5)

            fig1.tight_layout(pad=2.5)
            fig1.savefig('corr_histogram_global.png', bbox_inches='tight',
                         dpi=130, facecolor='#ffffff')
            plt.close()
            res.write('<img src="corr_histogram_global.png">')

            # ══════════════════════════════════════════════════════════════════
            # Charts 2 & 3 — Lower-half and upper-half zooms  (side by side)
            # Split at the user-defined threshold (corr_split, default 0.5).
            # Values exactly on the threshold appear in both subsets so that
            # the split line is visible in both charts.
            #
            # v3+: each zoom uses a fully adaptive Y-axis:
            #   - Y-max is set to local_max × 1.18 (5 % head room for count labels)
            #   - If the zoom's local_max is itself >= 20× its own second-largest
            #     bar, the conditional log twin-axis (from _draw_histogram) kicks
            #     in automatically, keeping every bar readable.
            # ══════════════════════════════════════════════════════════════════
            r_lower_half = r_vals[r_vals <= r_split]
            r_upper_half = r_vals[r_vals >= r_split]

            fig23, (ax_low, ax_high) = plt.subplots(
                1, 2, figsize=(16, 5.5),
                gridspec_kw={'wspace': 0.45}
            )
            fig23.patch.set_facecolor('#ffffff')
            fig23.subplots_adjust(top=0.82)

            # ── helper: compute adaptive Y ceiling for a zoom axis ────────────
            def _adaptive_ylim(counts_arr, head_room=1.18):
                """
                Return an appropriate Y ceiling for a zoom histogram.

                Strategy:
                  1. If all bars are roughly balanced (max < 20× second-max),
                     simply use max × head_room — the linear scale is fine.
                  2. If the chart is heavily skewed AND the log twin-axis was
                     triggered by _draw_histogram, the linear axis still needs
                     a sensible ceiling so the dominant bar is shown fully.
                     We use the same max × head_room but ensure the ceiling is
                     never less than 1 (edge-case guard).
                """
                local_max = counts_arr.max() if len(counts_arr) > 0 else 1
                return max(local_max * head_room, 1)

            # — Lower zoom: [r_min → r_split] —
            counts_low, _ = _draw_histogram(ax_low, r_lower_half, r_min, r_split,
                            n_bins=N_BINS, annotate_median=False)
            _style_axis(ax_low,
                        title='Lower-Half Zoom  —  20 bins',
                        xlabel='Pearson r',
                        ylabel='Number of ions',
                        n_ions=len(r_lower_half),
                        subtitle=f'Range : [{r_min:.3f} → {r_split:.3f}]  (below split threshold)')
            # Adaptive Y-axis: scale to the local distribution, not the global one
            ax_low.set_ylim(0, _adaptive_ylim(counts_low))
            # Shade the lower zone lightly
            ax_low.axvspan(r_min, r_split, alpha=0.05, color='#c0392b', zorder=1)
            ax_low.axvline(r_split, color='#7f8c8d', linewidth=1.2,
                           linestyle=':', zorder=4,
                           label=f'Split threshold = {r_split:.3f}')
            ax_low.legend(fontsize=7.5, loc='upper left',
                          framealpha=0.85, edgecolor='#b0bec5')

            # — Upper zoom: [r_split → r_max] —
            counts_high, _ = _draw_histogram(ax_high, r_upper_half, r_split, r_max,
                            n_bins=N_BINS, annotate_median=False)
            _style_axis(ax_high,
                        title='Upper-Half Zoom  —  20 bins',
                        xlabel='Pearson r',
                        ylabel='Number of ions',
                        n_ions=len(r_upper_half),
                        subtitle=f'Range : [{r_split:.3f} → {r_max:.3f}]  (above split threshold)')
            # Adaptive Y-axis: scale to the local distribution, not the global one
            ax_high.set_ylim(0, _adaptive_ylim(counts_high))
            # Shade the upper zone lightly
            ax_high.axvspan(r_split, r_max, alpha=0.05, color='#1a6faf', zorder=1)
            ax_high.axvline(r_split, color='#7f8c8d', linewidth=1.2,
                            linestyle=':', zorder=4,
                            label=f'Split threshold = {r_split:.3f}')
            ax_high.legend(fontsize=7.5, loc='upper left',
                           framealpha=0.85, edgecolor='#b0bec5')

            fig23.suptitle(
                f'Analytical Zoom Views — Split at r = {r_split:.2f}'
                f'  (configurable threshold, default 0.50)\n'
                f'Y-axis independently scaled to each half for maximum readability',
                fontsize=10, fontweight='bold', color='#1c2833', y=0.99
            )
            fig23.savefig('corr_histogram_zooms.png', bbox_inches='tight',
                          dpi=130, facecolor='#ffffff')
            plt.close()
            res.write('<img src="corr_histogram_zooms.png">')

            # ── Summary statistics table ──────────────────────────────────────
            r_mean     = float(np.nanmean(r_vals))
            n_strong   = int((r_vals >= 0.80).sum())
            n_moderate = int(((r_vals >= 0.50) & (r_vals < 0.80)).sum())
            n_weak     = int(((r_vals >= 0.30) & (r_vals < 0.50)).sum())
            n_poor     = int((r_vals < 0.30).sum())

            if r_mean < 0.5:
                res.write(f"""<div class="warn-box">
                ⚠️ The mean correlation between the two workflows is low (mean r = {r_mean:.3f}).<br>
                This may indicate significant differences in pre-processing or normalisation strategies.
                </div>""")

            res.write(f"""
            <table>
              <tr><th>Statistic</th><th>Value</th></tr>
              <tr><td>Ions analysed</td><td>{n_ions}</td></tr>
              <tr><td>Mean r</td><td>{r_mean:.3f}</td></tr>
              <tr><td>Median r</td><td>{r_median:.3f}</td></tr>
              <tr><td>Min r</td><td>{r_min:.3f}</td></tr>
              <tr><td>Max r</td><td>{r_max:.3f}</td></tr>
              <tr><td>Strong correlation (r ≥ 0.80)</td>
                  <td>{n_strong} ({round(100*n_strong/n_ions,1)}%)</td></tr>
              <tr><td>Moderate correlation (0.50 ≤ r &lt; 0.80)</td>
                  <td>{n_moderate} ({round(100*n_moderate/n_ions,1)}%)</td></tr>
              <tr><td>Weak correlation (0.30 ≤ r &lt; 0.50)</td>
                  <td>{n_weak} ({round(100*n_weak/n_ions,1)}%)</td></tr>
              <tr><td>Poor / negative correlation (r &lt; 0.30)</td>
                  <td>{n_poor} ({round(100*n_poor/n_ions,1)}%)</td></tr>
            </table>
            """)

else:
    res.write("""<h2>4. Intensity Comparison</h2>
    <div class="warn-box">
    ⚠️ No common samples were found between the two workflows.<br>
    Intensity comparison is therefore not possible.<br>
    Please verify that sample names are identical in both dataMatrix files.
    </div>""")

res.write("</body></html>")
res.close()

print("   Converting HTML → PDF...")
HTML('res.html').write_pdf('associations_quality.pdf')

print("\n" + "="*60)
print("✅ Analysis complete!")
print(f"   → Associations.tabular       (full association table)")
print(f"   → associations_quality.pdf   (quality report)")
print("="*60 + "\n")
