# MSOutputComparison

**Galaxy tool for comparing two mass spectrometry (MS) pipeline outputs**

Version: `1.5.0+galaxy0`

---

## 📋 Description

`MSOutputComparison` is a Galaxy tool that compares the results (ions, masses, retention times, intensities) of **two different MS pipelines**, potentially obtained from the same samples.
It addresses questions such as the following:

- Do variations in the pre-processing parameters influence the detected ions?
- Do my two sets of MS parameters yield equivalent results?
- What proportion of ions is shared between the two pipelines?

The tool uniquely couples ions from one pipeline to those of the other, based on a user-defined m/z and retention time (RT) tolerance,
then generates an association table and a quality report (in PDF format).

---

## 📦 Dependencies

| Package    | Version |
|------------|---------|
| numpy      | 2.1.0   |
| pandas     | 2.2.3   |
| matplotlib | 3.9.3   |
| seaborn    | 0.13.2  |
| weasyprint | 62.3    |

These dependencies are automatically managed by Galaxy via Conda/Bioconda (declared in the `<requirements>` section of the XML file).

---

## 📁 Repository files

- `MSOutputComparison.xml` — Galaxy wrapper (interface, parameters, help)
- `MSOutputComparison.py` — Python script performing the comparison and generating the report

---

## 🧬 Input files

For **each** MS pipeline, two files produced by the last step of the workflow (typically `groupChromPeaks` or `fillChromPeaks` under XCMS) are required:

1. **variableMetadata** — list of detected ions with their mass (m/z) and retention time (RT)
2. **dataMatrix** — intensities of each ion per sample

These tables should match the W4M's 3-tables format (see here for more information).

⚠️ To benefit from intensity comparisons, both pipelines must have been applied to at least a few **shared samples**.
Otherwise the tool will only produce the matching results, but without the intensity diagnosis.

---

## ⚙️ Parameters

| Parameter | Description | Default value |
|---|---|---|
| **Pipeline 1 / 2 name** | Short name used in the report (charts, titles) | `Pipeline 1` / `Pipeline 2` |
| **m/z tolerance** | Maximum mass distance (in Da) to match an ion from pipeline 1 with an ion from pipeline 2 | `0.02` |
| **RT tolerance** | Maximum retention time difference to match two ions (same unit as the `rt` column in the file, minutes or seconds) | `0.25` |
| **Correlation split threshold (r)** | Pearson r threshold used to split the two analytical zoom histograms in the report (section 5) | `0.5` |
| **m/z column name** | Exact name of the mass column in the `variableMetadata` file | `mz` |
| **RT column name** | Exact name of the retention time column in the `variableMetadata` file | `rt` |

### Example of tolerance values

- **m/z** — Orbitrap / QTOF (high resolution): 0.005–0.02 Da; Triple quadrupole (standard resolution): 0.05–0.5 Da
- **RT** — in minutes: 0.1–0.5; in seconds: 5–30

---

## 🧮 Matching algorithm

1. The user sets an m/z and RT tolerance.
2. For each ion in dataset 1, a list of candidates is built from the ions in dataset 2 located within the neighborhood defined by the tolerances.
3. For each pair (ion / candidate), the mass difference is calculated, producing a table `[ion dataset 1, candidate dataset 2, mass difference]`.
4. This table is sorted by ascending mass difference.
5. Two empty blacklists are initialized (one per dataset), along with an empty association table.
6. Starting from the smallest difference: if neither ion is already in its respective blacklist, they are paired and added to the blacklists; otherwise, move on to the next difference.
7. Ions left unmatched after this procedure are added to the table with the label `no match`.

This algorithm guarantees that each ion in a dataset is associated with the closest (by mass) ion in the other dataset.

---

## 📤 Output files

### `Associations.tabular`
Table listing all ions from both pipelines, including:
- ions shared by both pipelines (with their mass/RT differences)
- ions present only in pipeline 1 (column `idWF2 = no match`)
- ions present only in pipeline 2 (column `idWF1 = no match`)

The file is a text file with tabulation-separated values (TSV).

### `associations_quality.pdf`
Visual report containing:
- a matching summary (matched / unmatched ions)
- heatmaps of mass and RT differences (scale based on the tolerances entered)
- intensity comparison charts (when some samples are shared in both pipelines)
- three correlation histograms: global (20 bins, min→max), lower-half zoom (min→split threshold), upper-half zoom (split threshold→max), with automatic log axis in case of strong skew

---

## 💡 Usage tips

- **Too few ions matched** → gradually increase the m/z and RT tolerances to spot potentially threshold effects.
- **Almost all ions matched** with small differences → the two pipelines are very similar ✅.
- **Significantly different intensities** despite good matching → the pre-processing differs between the two pipelines.

---

## 🧪 Tests

One functional test is defined in the XML wrapper (`<tests>`), with default parameters (m/z = 0.02, RT = 0.25, corr_split = 0.5).

---

## 👤 Credits

- **Version 1.0.0**: Quentin Ruin, Mélanie Pétéra — INRAE / MetaboHUB
- **Version 1.5.0**: Abdou Lahat KA — INRAE

### What's new in version 1.5.0
- Correlation histograms with dynamic binning (step = (max−min)/20) and configurable lower/upper zooms
- Significantly improved visual design (professional color scheme, annotations)
- User-defined parameters for custom m/z and RT column names
- Removal of the strict naming requirement for the identifier column
- Reinforced input file validation (defensive programming)
- Heatmaps recalibrated on user-defined tolerances (no imposed bounds), with a fix for an m/z axis inversion bug
- Automatic dual Y-axis (linear + log) on histograms in case of strong skew between correlation classes

---

## 📚 Citation

DOI: [10.1093/bioinformatics/btu813](https://doi.org/10.1093/bioinformatics/btu813)

## 🆘 Support

France Bioinformatique community: https://community.france-bioinformatique.fr/c/galaxy/10
