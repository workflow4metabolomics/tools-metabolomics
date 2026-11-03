# PR #407 — ImputeMSmv: review consolidation and fix plan

Scope: [PR #407](https://github.com/workflow4metabolomics/tools-metabolomics/pull/407) — new tool
`tools/imputemsmv` (Galaxy wrapper for Delaporte *et al.* 2019 MAR/MNAR imputation).

Reviewed at commit `657d23acb` (== PR head).

Review inputs collected: 8 inline comments from `hechth` (CHANGES_REQUESTED), 1 from `lecorguille`,
9 inline comments from Copilot plus 2 additional items listed only in the Copilot review overview
(no inline comment). **No comment in the PR carries a `suggestion` block**, so every fix below is
proposed from scratch.

Every claim about code behaviour in this document was checked by executing the code in
`tools/imputemsmv` on R 4.3.3 against the PR's own test fixtures and against minimal synthetic inputs.
Where my verdict differs from a reviewer's, the measured behaviour is quoted. External facts (paper
identities, DOIs) were resolved against Crossref, not from memory — the two that appear are
`10.1016/j.chemolab.2019.03.005` and `10.18637/jss.v070.i01`.

---

## 0. Summary

| | Count |
|---|---|
| Reviewer comments triaged | 20 |
| — confirmed as described, fix needed | 16 |
| — real issue, but wrongly described (re-scope before fixing) | 3 |
| — refuted, no fix (reply and resolve) | 1 |
| **New findings from my own review** | **27** |
| — P0, corrupts scientific output silently | 4 |
| — P1, crash or wrong result on a reachable path | 4 |
| — P2, robustness / consistency / reproducibility | 13 |
| — P3, style, docs, lint | 6 |

Appendix A lists every reviewer comment with its verdict; Appendix B indexes all 27 new findings.

The two headline findings are not in any existing comment:

1. **P0-A — silent mislabelling of output columns and metadata rows.** The data matrix is rebuilt by
   concatenating per-group blocks, then its column names are overwritten *positionally*
   (`ImputeMSmv_function.R:333`). Whenever the replicate groups are not contiguous in the data
   matrix — which is the normal case for a randomised acquisition order — observed intensities are
   attached to the wrong samples. Variable metadata (`:336`) and sample metadata (`:339`) are
   `cbind`ed positionally too and misalign the same way. No test covers this; the fixtures use a
   perfectly blocked design.
2. **P0-B — the `min` MNAR method imputes from the wrong columns.** `min(Mat_noMNAR_temp[j, ])`
   (`:302`) takes its minimum over the bookkeeping columns `numero`/`NbNaN`/`pNbNaN`, not over
   intensities. On the PR's own test matrix, **86.6 % of imputed values land below 1** in a matrix
   whose observed intensities start at 402, and the largest imputed value is 2.1e8 (built from the
   feature row index). This is a user-selectable default-adjacent option with no test.

Both are consistent with the CI being green: 3 test blocks exercise exactly one parameter
combination (`lod` + `mean`), on blocked, VM-ordered input.

---

## 1. Responses to existing reviewer comments

Severity legend: **P0** corrupts scientific output silently · **P1** crash or wrong result on a
reachable path · **P2** robustness/consistency · **P3** style, docs, lint.

### 1.1 `hechth` — CHANGES_REQUESTED

---

#### C1 · `.shed.yml:2` — *"Please include some information on the method of imputation — either here or in the long description"*

**Status: accepted.** The `<description>` in the XML carries the method; the Tool Shed entry does not,
so the Tool Shed page is vaguer than the tool.

**Proposed fix** — extend `long_description` (and keep it consistent with the XML `<description>`):

```yaml
description: '[W4M][Preprocessing] Classifies and imputes MS missing values (MAR/MNAR)'
long_description: >
  Part of the W4M project: http://workflow4metabolomics.org.
  Classifies missing values per feature and per replicate group as MAR or MNAR from the
  within-group missingness rate, then imputes them with group-aware strategies: MAR from the
  group mean or median, MNAR from a global limit of detection or from a per-feature fraction
  of the lowest detected intensity, both combined with multiplicative uniform noise.
  Based on Delaporte et al. (2019), 10.1016/j.chemolab.2019.03.005.
```

Note: adding `citables: [true]` plus `citation` in `.shed.yml` is **not** the W4M convention here —
no other tool in this repo uses `citables`; citations live in the XML. Skip it.

---

#### C2 · `.shed.yml:5` — *"Can be removed as this should be mentioned explicitly in the citation section"*

**Status: accepted, with a caveat.** The intent is to stop describing the method as "based on
Delaporte et al. (2019)" in prose while the citations block points at the 2025 W4M platform paper.

That leaves a real gap, confirmed against Crossref:
`10.1016/j.chemolab.2019.03.005` = *Delaporte, Cladière, Camel, "Missing value imputation and data
cleaning in untargeted food chemical safety assessment by LC-HRMS", Chemom. Intell. Lab. Syst.* —
which **is** the method paper. So the method paper *is* already cited. Two consequences:

- Prose attributions to "Delaporte et al. (2019)" are correct, but the bibtex key is
  `delporte_workflow4metabolomics_2025` (the W4M platform paper) — the prose and the citation list
  are simply not linked. Use a `:cite:` role, or name the DOI in prose.
- The help text (XML `:246`) says "see Delaporte et al. (2019) for advanced arguments" while the
  `Contact:` line and the credits block credit **Cladière** as the original code author. Both are
  right (Cladière wrote the code, Delaporte *et al.* published the method), but the help never says
  which reference backs the *code*. Add one line distinguishing them.

**Proposed fix:** delete the "Based on the work presented by Delaporte et al. (2019)" sentence from
`.shed.yml` `long_description` (superseded by C1's wording, which names the DOI), and in the XML help
render the method reference as `:cite:`doi:10.1016/j.chemolab.2019.03.005`` so the prose resolves to
the bibtex/DOI entry.

---

#### C3 · `ImputeMSmv.xml:33` — *"This should be replaced by the detect_errors field on the command section"*

**Status: accepted.** `<stdio><exit_code range="1" level="fatal"/></stdio>` is legacy; 20+ tools in
this repo already use `<command detect_errors="exit_code">` (e.g. `w4mconcatenate/concat.xml:13`,
`batch_dispersion/batch_dispersion.xml`).

**Proposed fix:**

```xml
-  <stdio>
-    <exit_code range="1" level="fatal" />
-  </stdio>
-
-  <command>
+  <command detect_errors="exit_code">
```

---

#### C4 · `ImputeMSmv_wrapper.R:207` — *"rather use named lists for the result_tables — indices can often lead to confusion"*

**Status: accepted.** `imputemsmv()` already returns a named list
(`ImputeMSmv_function.R:341`: `list(Datamatrix=, sampleMetadata=, VariableMetadata=)`) — only the
consumer uses `[[1]]/[[2]]/[[3]]`.

**Proposed fix** in `ImputeMSmv_wrapper.R:200-207`:

```r
DM_out <- data.frame(FirstCol = rownames(d), result_tables$Datamatrix,      check.names = FALSE)
SM_out <- data.frame(FirstCol = rownames(s), result_tables$sampleMetadata,  check.names = FALSE)
VM_out <- data.frame(FirstCol = rownames(v), result_tables$VariableMetadata, check.names = FALSE)
colnames(DM_out)[1] <- dName
colnames(SM_out)[1] <- sName
colnames(VM_out)[1] <- vName
```

Two notes: (a) `"df"` as the placeholder column name is misleading (it is the identifier column) —
rename it and swap once, or build with `check.names=FALSE` and assign directly; (b) this comment also
asks for the change in `ImputeMSmv_function.R` — the returned list is already named, so the only
change needed there is the `@return` roxygen block (see R22).

---

#### C5 · `ImputeMSmv.xml:165` — *"What is actually being plot if none of the plots above is selected?"*

**Status: accepted — and the current answer is "an unusable 0-page PDF".** Verified behaviour:

- `plots` is `multiple="true"` with `optional="true"` and no `selected`; nothing chosen ⇒ Galaxy
  passes `--plots ''` ⇒ `strsplit("", ",")[[1]]` is `character(0)` ⇒ `"None" %in% chosen_plots` is
  `FALSE` ⇒ the `else` branch runs and logs *"plotting 0 chart(s)"*.
- `pdf(); dev.off()` with no plot emitted produces a **valid PDF with 0 pages** (`file` reports
  `PDF document, version 1.4, 0 page(s)`), which Galaxy presents as a green, empty dataset.
- The `"None"` branch at `ImputeMSmv_wrapper.R:223` is **dead code**: `None` is not among the option
  values in the XML (options are `DP HM BP1 BP2 PC`), so it is unreachable.

**Proposed fix** — make the intent explicit and always emit a renderable page:

```xml
<param name='plots' type='select' display='checkboxes' multiple='true' label='Plots to generate'
       optional='true' help='...If none is selected, the PDF will contain a single page stating that
       no plot was requested.'>
  <option value='None'>No plot</option>
  <option value='DP'>Density plot (before vs after imputation)</option>
  ...
```

```r
if (length(chosen_plots) == 0 || "None" %in% chosen_plots) {
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))   # plot.new() alone leaves no frame; be explicit
    title(main = "NO PLOT SELECTED")
}
```

Keep the guard `length(chosen_plots) == 0` — that is the reachable case. Add a test asserting the
PDF is non-empty (`compare="size"` with a minimum, or a checksum fixture) so this cannot regress.

---

#### C6 · `ImputeMSmv_function.R:88` — *"Rather do not supply defaults here, since the script is only used from galaxy… This could introduce errors later"*

**Status: accepted — and this is load-bearing, not just style.** `--mnar_computation` is emitted only
on the `lod` branch of the XML (`ImputeMSmv.xml:63-68`), so with **MNAR method = `min`** the argument
never arrives. Today `mnar_computation = c("quantile","mean","median")` silently resolves to
`"quantile"`; removing the default exposes it as an error. Measured:

```
f(mnar_computation = NULL)      -> "quantile"      # NULL also silently falls back to the default
f()  with default removed       -> ERROR: argument "mnar_computation" is missing, with no default
```

**Proposed fix** (do this together with R3 in §2 — the same three lines):

```r
imputemsmv <- function(DataMatrix, sampleMetadata, variableMetadata, colClass, LimNan,
                       setSeed, seed, mnar_method, percLOD, SDLOD, pmin, mar_method, SDmean,
                       mnar_computation, advanced, digits) {
    mnar_method      <- match.arg(mnar_method,      choices = c("lod", "min"))
    mar_method       <- match.arg(mar_method,       choices = c("mean", "median"))
    mnar_computation <- match.arg(mnar_computation, choices = c("quantile", "mean", "median"))
```

and pass `mnar_computation` unconditionally from the command block, giving the `min` branch a
harmless value (or better: drop the parameter entirely from the `min` path — see §2 R3).

---

#### C7 · `ImputeMSmv_function.R:101` — *"This actually won't work, since the sample is based on the seed… generating a random number to use as a seed does not work"*

**Status: accepted in substance; one mechanism detail corrected.** `sample(1:10000, 1)` *is*
correctly drawn from a fresh RNG state (R seeds automatically), so the printed value is a real,
reusable seed — a user who types it back **does** reproduce the run *of the seed draw*, but not of
the imputation, because nothing seeds the RNG before the imputation draws consume it. The comment's
conclusion is right: the printed number is not actionable in the `setSeed = FALSE` case. Measured
additionally: `set.seed("123")` and `set.seed(123)` produce identical streams, so the string seed is
harmless (see C-10 for where that matters instead).

**Proposed fix** — either make the printed value honest, or drop the claim:

```r
if (isTRUE(setSeed)) {
    set.seed(seed)
    cat("Seed set by the user:", seed, "\n")
} else {
    # Seed from the clock so the run is non-reproducible by design, and report what we actually did.
    random_seed <- as.integer(Sys.time() * 1000) %% .Machine$integer.max
    set.seed(random_seed)
    cat("No user seed. RNG seeded from the clock:", random_seed,
        " - record it to reproduce this exact run.\n")
}
```

Then **delete the duplicate `set.seed()` in the wrapper** (`ImputeMSmv_wrapper.R:180`) — seeding
twice, in two files, for one option is how the two drift apart. Keep seeding in exactly one place
(the function, which owns the draws). Also update the XML `help` for `repro.enable`, which currently
promises "the seed used will be outputed in the log and thus can be used for reproducibility".

---

#### C8 · `ImputeMSmv_function.R:110` — *"This validation should be part of the wrapper, not the function… only run the computations if everything is ok rather than starting computations and then checking"*

**Status: accepted.** Note the line has drifted: the empty-feature check is now at `:108-110`. It sits
after 100 lines of logging and inside the group loop that does the actual work, so it can only fire
after the matrix has been split by group.

**Proposed fix** — move every precondition to a single wrapper-side validator that runs immediately
after import, before `imputemsmv()` is called. This is the natural home for the checks that are
currently scattered or missing, i.e. §2 findings R1–R7 and C-6. Proposed shape (in
`ImputeMSmv_wrapper.R`, replacing the ad-hoc blocks at `:85-145`):

```r
validate_inputs <- function(d, s, v, col) {
    problems <- character(0)
    require_ <- function(cond, msg) if (!isTRUE(cond)) problems <<- c(problems, msg)

    require_(nrow(d) > 0,    "The data matrix contains no features.")
    require_(ncol(d) > 0,    "The data matrix contains no samples.")
    require_(nrow(v) == nrow(d), "The data matrix and the variable metadata describe a different number of features.")
    require_(nrow(s) == ncol(d), "The data matrix and the sample metadata describe a different number of samples.")
    require_(setequal(rownames(v), rownames(d)), "Feature identifiers differ between the data matrix and the variable metadata.")
    require_(setequal(rownames(s), colnames(d)), "Sample identifiers differ between the data matrix and the sample metadata.")
    require_(all(vapply(d, is.numeric, logical(1))), "The data matrix contains non-numeric values.")
    require_(col %in% colnames(s), paste0("Column '", col, "' is not found in the sample metadata."))
    require_(!any(is.na(s[[col]])), paste0("Column '", col, "' contains missing values."))
    require_(sum(is.na(d)) > 0, "No missing values found in the data matrix; nothing to impute.")
    if (length(problems) == 0) {
        fully_missing <- rownames(d)[rowSums(is.na(d)) == ncol(d)]
        require_(length(fully_missing) == 0,
                 paste0("These features are missing in all samples, so no reference value exists: ",
                        paste(fully_missing, collapse = ", ")))
    }
    if (length(problems)) stop(paste(problems, collapse = "\n"), call. = FALSE)
}
```

and reduce the function to imputation plus one defensive `stopifnot`. Also reject the
`colClass` values that break the object-name scheme (§2 R6) here, where the error message can name
the offending characters.

---

#### C-4b · `ImputeMSmv_function.R:82,86` — same class of issue as C6

`mnar_method = c("lod","min")` and `mar_method = c("mean","median")` are defaults too, but they are
harmless *because* the XML always passes both. Still, drop them for consistency with C6 and pass
`choices=` explicitly, so a missing parameter is always an error rather than a silent first option.

---

### 1.2 `lecorguille`

---

#### C9 · `ImputeMSmv.xml:219` — *"Could it be replace by the `creator` section?"* (the `**Contact:**` block in the help)

**Status: accepted.** The help hard-codes contact details that duplicate — and disagree with — the
`<creator>` block at `:15-19`: help says "Mélanie Pétéra", creator says `givenName='Anthony'`;
`<creator>` carries no email at all, so it cannot serve as the contact route. 10 other tools in this
repo use `<creator>`, several with `email=`.

**Proposed fix** — delete the `**Contact:**` paragraph (`ImputeMSmv.xml:217-219`) and route support
through the machine-readable block plus the help's existing support URL:

```xml
<creator>
  <person givenName="Anthony"   familyName="Sleiman"  jobTitle="Internship student" />
  <person givenName="Mathieu"   familyName="Cladiere" jobTitle="Researcher" />
  <person givenName="Melanie"   familyName="Petera"   jobTitle="Engineer" email="pfem-tools@inrae.fr" />
  <organization name="PFEM - INRAE - MetaboHUB" url="https://eng-pfem.isc.inrae.fr/" />
  <organization name="Workflow4Metabolomics"    url="https://workflow4metabolomics.org/" />
</creator>
```

Use ASCII in `<creator>` attributes (the help text uses "Pétéra" and "Cladière"); keep the accented
forms only in the prose credits. Confirm the actual support address with the maintainer before
filling `email=` — do not invent one. Keep the existing
`https://community.france-bioinformatique.fr/c/workflow4metabolomics/10` link (currently only in
README) in the help as the user-facing support route.

---

### 1.3 Copilot

---

#### C10 · `ImputeMSmv_function.R:226` — *"MAR rows are selected only by `pNbNaN <= LimNan`… one missing value in a two-replicate group at `LimNan = 0.5` is classified MNAR yet is included in both this MAR table and the MNAR table, duplicating features during the `rbind` below. Select rows from the computed `MARMNAR` labels instead."*

**Status: real bug, wrong mechanism and wrong example — re-scope before fixing.** Measured:

- The two selections are **complementary**, not overlapping: MAR is `pNbNaN <= LimNan` (`:226`), MNAR
  is `pNbNaN > LimNan` (`:278`, `:297`). No row can be in both. On the real test matrix the
  `rbind` at `:321` returns exactly `nrow(DataMatrix)` rows with no duplicates.
- The stated example cannot even be produced by this code: a group of 2 replicates has
  `rep_num = 2 < repTH`, so `MARMNAR` is forced binary and `LimNan` is not consulted at all.
- The **actual** defect is the mirror image of the one described, and is worse than a duplicate:
  the *label* comes from `MARMNAR` (which applies the small-group binary rule) while the *strategy*
  comes from a raw `pNbNaN` comparison (which does not). So a feature labelled **MNAR** in the output
  variable metadata is silently imputed with the **MAR** strategy. Reproduced:

  ```
  group g, 2 replicates, LimNan = 0.5, mnar_method = "min", pmin = 0.8
  feature fMAR: 1 of 2 values missing  ->  pNbNaN = 0.5
    reported class : MNAR          (binary small-group rule)
    imputed value  : 1000          (= group mean, i.e. the MAR strategy)
    expected       : ~0.8 * 1000 = 800 (the MNAR strategy the user selected)
  ```

  This is reachable **without** the advanced option: `LimNan1` defaults to 0.25 but allows up to 0.5,
  and any group with 2 replicates and `pNbNaN <= LimNan` hits it.

**Proposed fix.** Use the computed label as the single source of truth, in all three places. Replace
`:165-190` with a vectorised form and select on it:

```r
MARMNAR <- ifelse(pNbNaN == 0, "No_MV",
             ifelse(rep_num >= repTH & pNbNaN <= LimNan, "MAR", "MNAR"))
...
mar_rows  <- MARMNAR == "MAR"
mnar_rows <- MARMNAR == "MNAR"
assign(nom_DataMatrix_gp_MAR[i],  gp_mat[mar_rows,  , drop = FALSE])
assign(nom_DataMatrix_gp_MNAR[i], gp_mat[mnar_rows, , drop = FALSE])
```

and drop the `subset(..., pNbNaN <= LimNan)` / `subset(..., pNbNaN > LimNan)` calls entirely. The
labels and the strategies then cannot disagree by construction, and the `No_MV` features (currently
dropped from *both* subsets and re-introduced only because the MAR/MNAR blocks happen to keep their
input values) should be passed through explicitly — see R10 in §2.

---

#### C11 · `ImputeMSmv_function.R:259` — *"`DataMatrix` is a data frame here… base R's default quantile method requires a numeric vector, so the default `lod` mode—and the XML test—fails before any imputation. Flatten the matrix before calculating the quantile."*

**Status: REFUTED. Do not change; reply and resolve.** `quantile.default` calls `stat.format <-
stats:::.calc.ucrt` on `unlist`-compatible input and flattens data frames correctly. Measured:

```r
d <- data.frame(a=c(1,2,NA,4), b=c(100,200,300,400))
length(quantile(d, na.rm=TRUE, prob=0.03, names=FALSE))   # 1        <- flattened, not column-wise
quantile(d, na.rm=TRUE, prob=0.03, names=FALSE)           # 1.18     == quantile(unlist(d), ...)
# column-wise would have been 1.06 and 109
```

Registered methods are `default`, `ecdf`, `POSIXt` — a `data.frame` dispatches to `default`, which
handles it. And empirically the default `lod` mode runs end-to-end on the PR fixtures and
reproduces `outDmShT1.tabular` byte-for-byte (§3.1). The suggested `unlist()` is a cosmetic no-op;
adding it would imply a defect that does not exist. **No action.**

---

#### C12 · `ImputeMSmv_function.R:302` — *"`Mat_noMNAR_temp` still contains the bookkeeping columns `numero`, `NbNaN`, `pNbNaN`, so this minimum is often taken from a row index or proportion rather than an intensity. It also only examines the current group, whereas the documented method uses the per-feature minimum across all samples; use the corresponding original row in `DataMatrix`."*

**Status: confirmed on both counts — this is P0-B. Highest-priority Copilot finding.** Measured on the
PR's own fixtures with `mnar_method="min"`:

```
observed intensity range  : [402.1, 9.27e+08]
imputed (min) range       : [0.370, 2.09e+08]
imputed values < 1        : 6346 of 7325  (86.6 %)   — all 6346 below the observed minimum
imputed values < 1e-3 * min(o): 8
for contrast, lod imputed range: [~719, ~1123]
```

and in isolation, with the noise forced to exactly 1 (`SDLOD=0`, `pmin=0.8`):

```
feature f3, 4 replicates, observed 3000 3100 -> min over MNAR row incl. bookkeeping cols = 0.5
imputed f3[s3] = 0.4     expected 0.8 * 3000 = 2400
feature f2, imputed f2[s2] = 2036.667   expected 0.8 * 2000 = 1600
```

The reference for one feature is `min(numero, intensities…, NbNaN, pNbNaN)`, i.e. `min(1, …)` for
features with row index 1 and `min(pNbNaN, …)` whenever the proportions are smaller than the
intensities — and the 2.1e8 outliers come from `numero` for high-index features. The second half of
the comment is also right: the documented method ("the minimum of non-missing values is computed"
per feature, `ImputeMSmv.xml:306-307`, and the wrapper's own log at `:172`, "the lowest value
computed per feature **across all samples**") uses all samples, whereas `:302` sees only the current
group's columns.

**Proposed fix** — compute the reference from the untouched matrix, over all samples, once per
feature:

```r
# per-feature minimum over ALL samples, computed before any bookkeeping column exists
feat_min <- apply(as.matrix(DataMatrix), 1, function(x) min(x, na.rm = TRUE))
...
for (i in seq_along(groups)) {
    sub <- gp_mat[mnar_rows, , drop = FALSE]
    ref <- pmin * feat_min[rownames(sub)]
    ...
}
```

Guard the degenerate cases found in §2 R1/R4: `min(numeric(0), na.rm=TRUE)` returns `Inf` with a
warning, so a feature with no observed value anywhere must be rejected by the wrapper validator
(C8) *before* reaching here, and `feat_min` must never be allowed to be `Inf` or `NA`.

---

#### C13 · `ImputeMSmv_wrapper.R:221` — *"`on.exit()` is called at top level, outside any enclosing R function, so this aborts the wrapper before plotting; even inside a function it would register a second `dev.off()` after the explicit close at line 261. Remove this registration and close the PDF exactly once."*

**Status: wrong on both mechanisms, right conclusion — fix it, but as a lint/hygiene item, not a
crash.** Measured on R 4.3.3:

```
pdf(file=...); on.exit(dev.off(), add=TRUE); cat("reached\n"); dev.off()
-> "reached", exit code 0, no error, no warning
```

It does not abort, and it does not double-close: at top level the handler is registered into
`R_GlobalContext`, but `Rscript` exits via a longjmp path that never runs it (verified by making the
handler write a marker file — no marker was created). The `dev.off()` at `:261` is the only close
that executes. So the severity is **P3**, not critical.

**Proposed fix** — remove the line as asked, and wrap plotting in a function so the cleanup actually
works and any `stop()` inside a plot still closes the device (this matters: a plotting error today
leaves Galaxy with a truncated PDF and a confusing secondary error):

```r
render_plots <- function(path, width = 11, height = 7, expr) {
    con <- pdf(file = path, width = width, height = height)
    on.exit(invisible(dev.off(con = con)), add = TRUE)   # inside a function: this does run
    force(expr)
}
```

Note for the reply: also flag that `dev.off()` at `:261` writes `null device 1` into the job stdout
for every run; `invisible()` or `graphics.off()` avoids the noise.

---

#### C14 · `ImputeMSmv_function.R:102` — *"When `setSeed` is false, this value is only printed; it is not used to seed the RNG… the advertised logged seed cannot reproduce the run. Call `set.seed(random_seed)` before the imputation draws (or log a reproducible RNG state instead)."*

**Status: confirmed — same defect as C7, from the other side.** This is the accurate description of
the mechanism; C7's "does not work" reasoning is weaker but its conclusion matches. Treat C7 and C14
as one fix (see C7). Severity **P2** — reproducibility contract, not numeric corruption, and the
default XML help already warns results "will vary slightly".

One extra wrinkle the comment misses and the fix must handle: the wrapper *also* prints a claim about
a random seed (`ImputeMSmv_wrapper.R:182`) while doing nothing, so both places must be reconciled
with whichever behaviour is chosen. Prefer the clock-seed variant in C7: it keeps the documented
feature ("record the seed to reproduce") true.

---

#### C15 · `ImputeMSmv_function.R:590` — *"This normalization divides by the total count of each category without handling a zero total. A valid dataset with no MAR (or no MNAR/No_MV) entries therefore produces `NaN` values."*

**Status: confirmed, and worse than described — the plot draws literal `NaN%` labels.** Measured:

```r
No_MV <- c(5,5); MAR <- c(0,0); MNAR <- c(1,2)
cbind(No_MV = No_MV/sum(No_MV)*100, MAR = MAR/sum(MAR)*100, MNAR = MNAR/sum(MNAR)*100)
#      No_MV MAR     MNAR
# [1,]    50 NaN 33.33333
# [2,]    50 NaN 66.66667
# labels drawn by the code: "50%" "50%" "NaN%" "NaN%" "NaN%" "NaN%"
```

`barplot()` tolerates the `NaN` column silently, so the user sees a barplot annotated `NaN%` — no
error, no warning. Reachable in the default configuration: any dataset where one category is absent
overall. (This is not the same as "BP2 plot fails" as the overview claims — it does not fail; it
renders wrong annotations, which is more dangerous.)

**Proposed fix** — a small shared helper (also used by `plt_group`, `plt_pie`, §2 R11):

```r
# proportions of `x` as percentages; a zero total yields all-zero rather than NaN
pc <- function(x) { s <- sum(x); if (!is.finite(s) || s == 0) rep_len(0, length(x)) else x / s * 100 }
```

plus: skip drawing a category column whose total is zero, and label the axis accordingly, rather
than showing an empty stack with `%` labels.

---

#### C16 · `ImputeMSmv_function.R:661` — *"The filtering mask is recomputed after `info` has already been shortened. When a category is absent, `color` and `label` are then filtered using positions from the shortened vector, so the remaining slices receive incorrect colors/legend labels (for example, MNAR can be labeled MAR). Save the three-element mask before subsetting any vector."*

**Status: confirmed — real bug, but the observable symptom differs from the one described, so the
fix must be verified against the real symptom.** The mechanism is exactly as stated (`:659-661`:
`info <- info[info > 0]` then `color[info > 0]`). The detail that changes the outcome: R **recycles**
a logical subscript shorter than the target, so `color`/`label` come back at full length rather than
truncated:

```r
c("a","b","c")[c(TRUE, TRUE)]        # "a" "b" "c"   <- recycled, not 2 elements

true composition: No_MV 50 %  MAR 0 %  MNAR 50 %
  slices passed to pie() : 2   (50, 50)
  colours passed to pie(): 3   (white, PURPLE, MAROON)  -> slice colours actually used: white, PURPLE
  legend entries         : 3   (No_MV, MAR, MNAR)
```

So the real, user-visible defects when any category is absent are: (i) the **legend lists categories
with 0 %**, (ii) slices get the **wrong colour** (the MNAR slice is drawn in the MAR colour, because
`pie()` takes the first *n* colours), and (iii) with a single non-empty category the legend still
shows three entries. "MNAR labelled MAR" is not what happens — the labels stay aligned, the colours
do not.

**Proposed fix** — as asked, take the mask once; and use a fixed category→colour map so alignment
cannot depend on vector length:

```r
CATS  <- c(No_MV = "white", MAR = "#6551CC", MNAR = "#B22C2C")
counts <- c(No_MV = sum(groups == "No_MV"), MAR = sum(groups == "MAR"), MNAR = sum(groups == "MNAR"))
keep   <- counts > 0                                  # single mask, computed once
info   <- pc(counts[keep]); col <- unname(CATS[keep]); lab <- names(counts)[keep]
pie(info, labels = paste0(round(info, 2), "%"), col = col)
legend("topright", legend = lab, fill = col, border = "black", cex = 0.6)
```

`plt_group` (`:526-529`) and `plt_type` (`:596`, `:610`) should use the same fixed map. Checked
`plt_type` specifically for the same class of error: it colours the bars with
`hcl.colors(nrow(info))` and the legend with `hcl.colors(length(groups))`, which looked like two
different quantities. Measured, both are the number of groups (`info` is groups × categories, and
`groups` is the data frame of one column per group), so bar colours and legend colours do agree —
**not a bug**. The objection to `plt_type` is only that the same number is derived twice by different
routes; use one variable and the shared map.

---

#### C17 · `ImputeMSmv_wrapper.R:115` — *"This check operates on `sm_new`, which is a single row of per-sample NA counts, so it cannot detect a feature row that is entirely missing… test the row-wise NA counts of `d` instead."*

**Status: confirmed, and it is silently broken on the PR's own test data — not merely "can produce
invalid Inf references".** Measured on `inDmShT1.tabular`:

```
features that are 100 % NA in the fixture : M285T579, M331T579, M332T579, M353T579   (4 features)
sm_new  <- rbind(apply(d, 2, function(x) sum(is.na(x))))
rownames(sm_new)                          : NULL          <- the single row is unnamed
row_na  <- rownames(sm_new)[apply(sm_new == 0, 1, all)] : NULL
=> guard fires? FALSE
```

Two independent defects: it looks in the wrong direction (samples, not features), and `rbind` of a
`1×n` named vector yields one unnamed row, so `rownames(sm_new)` is `NULL` and the guard is
**inert for every input**. The correct check is trivially available and does fire:
`sum(rowSums(is.na(d)) == ncol(d))` → 4.

The consequence today: those 4 features sail through, and in the `lod` path they are imputed from the
global LOD (defensible). In the `min` path, combined with C12, `min(..., na.rm=TRUE)` over an
all-NA row returns `Inf` with a warning, and `round(pmin * (Inf * 1), 2)` = `Inf` is written to the
output matrix. `min` on the fixtures happens not to produce `Inf` only because the reference is being
taken from `numero` (C12) — so **fixing C12 without fixing C17 will start emitting `Inf`s**. Fix both
together.

**Proposed fix:** fold into the validator from C8 (`rowSums(is.na(d)) == ncol(d)`), reporting the
offending feature names, and delete `ImputeMSmv_wrapper.R:114-118`.

---

#### C18 · `README.md:8` — *"The documented Galaxy ID is `imoutemsmv`, but the tool XML declares `id="imputemsmv"`."*

**Status: confirmed. Trivial.** **Proposed fix:** `* **@galaxyID**: imputemsmv`. While in the file, the
`@version` line says `1.0.0` but the XML says `1.0.0+galaxy0` — align them. Two more README items
while in the file: `@main usage` duplicates the XML `<description>` (keep one authoritative wording,
ideally by referencing the tool page), and the *Test files* provenance note ("derived from publicly
archived datasets available on SoftwareHeritage", `README.md:48-49`) gives no accession or archive URL
— link the same Software Heritage origin the help already cites (`ImputeMSmv.xml:246`) so the fixtures
can be re-derived.

---

#### C19 · Copilot overview only (no inline comment): *"use `para$mnar_methods` for logging"*

**Status: confirmed — and the reason it looks fine is a latent trap.** `ImputeMSmv_wrapper.R:164`
reads `para$mnar_method`, but the flag is `--mnar_methods` (`ImputeMSmv.xml:62`). It resolves anyway
because `$` on a list does **partial matching** when the exact name is absent:

```r
p <- list(mnar_methods = "lod"); p$mnar_method      # "lod"   <- partial match, no error
```

This is fragile: it breaks the moment any *other* key starting with `mnar_method` is added (e.g.
`--mnar_method_x`), and it silently returns `NULL` in that case, which then throws
"argument is of length zero" inside `if`. **Proposed fix:** use `para$mnar_methods` at `:164` and
`:230`, and adopt a single spelling for the XML flag, the R key and the function argument (see R9 —
the same drift exists as `mar_methods`/`mar_method`).

---

#### C20 · Copilot overview only: *"set `quote = FALSE` for the data-matrix output"*

**Status: confirmed as an inconsistency; P3.** `ImputeMSmv_wrapper.R:210-212` writes `SM_out`/`VM_out`
with `quote = FALSE` and `DM_out` with the default `quote = TRUE`, so feature identifiers come out as
`"M60T53"` — visible in the committed fixture `outDmShT1.tabular`, whose first line is
`"name"\t"Condition1_R1"…`. `w4mconcatenate/concatenation_wrapper.R:38-40` uses `quote = FALSE` for
all three, so the repo convention supports the change.

**Proposed fix:** add `quote = FALSE` at `:210`. Caveat to state in the reply: with `quote = FALSE` an
identifier containing a tab would corrupt the file — that is what `fix_head()`'s duplicate/NA checks
partly guard, and W4M identifiers do not contain tabs. Regenerate the three fixtures after this
change (see §3.1 — they are reproducible, so this is safe).

---

## 2. New findings from my own review

Ordered by severity. Each is measured, not inferred.

### P0 — silently corrupts scientific output

---

#### R1 · Output columns are mislabelled whenever replicate groups are not contiguous

`ImputeMSmv_function.R:332-333`:

```r
DataMatrix_NoMV <- cbind(as.data.frame(mget(nom_DataMatrix_gp_noMV[seq_along(nom_DataMatrix_gp)])))
colnames(DataMatrix_NoMV) <- colnames(DataMatrix)
```

Each group block carries its own sample columns; `cbind` concatenates them **in group order**, then
`colnames()` overwrites the names **positionally**. If the groups are contiguous in the data matrix
columns, group order == column order and the result is correct. Otherwise the names are shifted.
Reproduced on a W4M-valid input where the wrapper's own order check *passes*, so no sorting saves it:

```
data matrix columns : S1 S2 S3 S4          sample metadata rows: S1 S2 S3 S4  (identical => no sorting)
groups              : gA = {S1,S3}, gB = {S2,S4}      <- interleaved (randomised acquisition)
concatenated order  : S1 S3 S2 S4

input  f1 :   1  NA   3  NA
output f1 :   1   3 0.8 0.8     <- 3 (a genuine S1/S3 observation) is now reported under S2
expected   :   1 0.8   3 0.8

observed (non-NA) cells preserved:  CURRENT = FALSE,  reindexed = TRUE
```

The wrapper's guard at `:131-138` only compares `identical(rownames(s), colnames(d))`; when it fails
it sorts both **by name** — which does *not* make group order match column order, and when it passes
it does nothing at all. Randomised/instrument-order acquisition interleaves conditions constantly, so
this is the normal case, not an edge case. The blocked design of the fixtures (all `Cond1_*` columns
contiguous) is exactly why CI is green.

**Fix (minimal, and validated):** reindex by name instead of renaming positionally, at `:333`:

```r
DataMatrix_NoMV <- DataMatrix_NoMV[, colnames(DataMatrix), drop = FALSE]
colnames(DataMatrix_NoMV) <- colnames(DataMatrix)
```

Verified: `observed cells preserved` flips from `FALSE` to `TRUE` on the reproduction above. The
`drop = FALSE` matters (§2 R2). The robust version sorts groups deterministically and asserts
`setequal(colnames(DataMatrix_NoMV), colnames(DataMatrix))` before the assignment, so a silent
mismatch becomes a hard error.

---

#### R2 · Variable metadata and sample metadata are `cbind`ed positionally

Same class as R1, two more places:

- `:336` `variableMetadata_MV <- cbind(variableMetadata, VariableData_MARMNAR)` —
  `VariableData_MARMNAR` (`:213`) is built in **data-matrix row order**; `variableMetadata` is in
  whatever order the user's file has. `imputemsmv()` never validates or reindexes against VM row
  names. Reproduced with VM rows reversed:

  ```
  ground truth (per feature, per group MV count > 0 ?):  fX: gA=0 gB=1   fZ: gA=1 gB=0
  VM in DM order      ->  fX No_MV/MNAR   fY No_MV/No_MV   fZ MNAR/No_MV     (correct)
  VM rows reversed    ->  fZ No_MV/MNAR   fY No_MV/No_MV   fX MNAR/No_MV     (each row carries
                                                                             another feature's labels)
  ```
  The output is internally consistent (row names come from VM), so nothing downstream notices: the
  user's per-feature MV classification is silently assigned to the wrong features.
- `:339` `sampleMetadata_MV <- cbind(sampleMetadata, SM_new)` — `SM_new` (`:128`) is computed from
  `apply(DataMatrix, 2, …)`, i.e. in **data-matrix column order**, and attached to sample metadata by
  position. With `rownames(s) = S3 S4 S1 S2` and true per-sample MV counts `S1=2, S2=2, S3=0, S4=0`,
  S3/S4 are reported as `66.67`.

**Fix:** key everything by name and never rely on positional `cbind` across the three tables:

```r
VariableData_MARMNAR <- VariableData_MARMNAR[rownames(variableMetadata), , drop = FALSE]
variableMetadata_MV  <- cbind(variableMetadata, VariableData_MARMNAR)

SM_new <- SM_new[rownames(sampleMetadata), , drop = FALSE]
sampleMetadata_MV <- cbind(sampleMetadata, SM_new)
```

plus a `setequal()` assertion. Better still: name `SM_new` when computing it
(`setNames(..., colnames(DataMatrix))`) so the alignment is checked rather than assumed. This belongs
with the C8 validator, which should reject differing *sets* and let the function reindex differing
*orders*.

---

#### R3 · The `min` MNAR method is broken (C12) — see C12 for the fix

Listed here so it is not lost: a user-selectable method produces imputed values two to four orders of
magnitude below the data's dynamic range on 86.6 % of cells, with no test. Add a test that asserts the
imputed values fall inside `[pmin * min(all observed), pmin * max(all observed) * (1 + SDLOD)]`.

---

#### R4 · The all-NA-feature guard is inert (C17) — see C17 for the fix

Inert for every input (`rownames(sm_new)` is `NULL`), and the fixtures contain 4 such features.
Interacts with R3/C12: after fixing C12 this becomes an `Inf` producer.

### P1 — crashes on reachable paths

---

#### R5 · Crash when any group has exactly one MAR feature

`ImputeMSmv_function.R:227-228`:

```r
Mat_noMAR_temp <- t(get(nom_DataMatrix_gp_MAR[i]))
Mat_noMar_temp_2 <- Mat_noMAR_temp[!(row.names(Mat_noMAR_temp) %in% c("NbNaN", "pNbNaN")), ]
```

The MAR subset of a group has columns `numero`, one per sample in the group, `NbNaN`, `pNbNaN` — i.e.
`n_reps + 3` columns. With exactly **one** MAR feature that subset has one *row*, so `t()` produces an
`(n_reps + 3) × 1` object, and the row-name filter at `:228` indexes it with a single column left over —
which R drops to a plain numeric vector (verified on a 2-replicate group: `class = numeric`,
`dim = NULL`, `nrow = NULL`, three elements left). So the trigger is the single *column*, not the row
count. `:232`/`:238` then evaluate `seq_len(nrow(...))` = `seq_len(NULL)`, which warns *"first element
used of 'length.out' argument"* and errors *"argument must be coercible to non-negative integer"*.

Reproduced on the default configuration, no advanced option needed.

*Corrected after re-measurement while writing the regression tests.* The exact trigger is:
**a group whose MAR subset contains exactly one feature.** Two of the three examples originally
listed here did not actually reproduce it, and why they don't is itself informative:

```
1 feature, 1 group of 4 replicates, 2 MVs (pNbNaN = 0.5 <= LimNan 0.5), mar_method = mean
  -> CRASH: Error in seq_len(nrow(Mat_noMar_temp_2)):
            argument must be coercible to non-negative integer

no-replicate design (3 samples, 3 distinct group labels), 2 features, where exactly
one feature is No_MV in the singleton group -- both mnar methods
  -> CRASH, same error

3 samples / 3 groups where BOTH features are NA in the same sample
  -> NO crash: the MAR subset is empty, nothing is indexed

2 groups of 2 replicates, one all-NA feature (pNbNaN = 1.0) and one complete feature
  -> NO crash: the all-NA feature is *labelled* MNAR by the small-group rule, but the
     row SELECTOR at :226 still uses `pNbNaN <= LimNan`, so it lands in the MAR subset
     anyway, giving two rows. R10/C10 masks R5 here.
```

That last point matters for the fix order: because `:226` selects by `pNbNaN <= LimNan`
rather than by the computed `MARMNAR` label, features the tool *calls* MNAR are swept through
the MAR loop (R10). Once C10/R18 make the selection honour the label (step 3), which subsets
have exactly one row changes, so R5's reachable cases shift. `drop = FALSE` is required either
way, and it is what makes step 3 safe to attempt.

The no-replicate case is the most serious: **the single-sample-per-group design documented in
the tool's own help text crashes the tool** ("If the design has no replicates, the column
should contain distinct values for all samples", `ImputeMSmv.xml:91`).

**Fix — one word, and it is behaviour-preserving.** Add `drop = FALSE` to the row filter at `:228`:

```r
Mat_noMar_temp_2 <- Mat_noMAR_temp[!(row.names(Mat_noMAR_temp) %in% c("NbNaN", "pNbNaN")), , drop = FALSE]
```

Verified on the two-feature case that previously worked: same seed, `all.equal` result is `TRUE`, so
nothing that runs today changes; and the one-feature case stops collapsing to a vector
(`CRASH: collapsed to a vector` → `ok`). Apply the same to the `t()` at `:227` and to the MNAR-side
subsets for symmetry.

Do **not** "tidy" this by dropping the transpose and iterating features outer — measured, that flips
the RNG consumption order and changes *every* imputed value for a given seed (10 of 10 NA cells differ
on a 4×6 block; see R21). Keep the transpose until step 8, where the value change is deliberate,
announced, and accompanied by regenerated fixtures.

Add a regression test with exactly one MAR feature in a group (T3).

---

#### R6 · `mnar_computation` never reaches the function on the `min` path

Consequence of C6/C-4b, promoted to a separate item because removing the defaults (as C6 asks) turns
this from silent into fatal *today* — i.e. **the fix for C6 breaks the `min` method unless the command
block is fixed in the same commit.** `ImputeMSmv.xml:63-68` emits `--mnar_computation` only under
`mnar_method == 'lod'`; the wrapper passes `para$mnar_computation` unconditionally (`:192`), which is
`NULL` for `min`.

**Fix (choose one, in the same commit as C6):** (a) always emit `--mnar_computation` by moving it out
of the conditional; or (b) preferred — drop the parameter from `imputemsmv()` when
`mnar_method == "min"`, since it is meaningless there, and guard the `match.arg` accordingly. Note
`mnar_computation` is currently used *only* inside the `lod` branch (`:260-269`), so (b) is a clean
removal.

---

#### R7 · `percLOD`/`pmin` share the same problem

Same XML conditional shape (`:63-68`) puts `--percLOD` on the `lod` branch and `--pmin` on the `min`
branch, while the wrapper passes both unconditionally (`:191`). Today they are used only inside their
respective branch, so `NULL` is harmless — but it means the function signature can receive `NULL` for
a documented numeric parameter, which will blow up the moment any default-free validation is added.
Resolve all three parameters the same way as R6.

#### R8 · `LimNan2` has `value=''` and no `required`, silently collapsing the whole matrix to MNAR

`ImputeMSmv.xml:154`:

```xml
<param name='LimNan2' type='float' value='' min='0.50' max='1' label='Extended missing value rate threshold' ... />
```

Two problems: the parameter is **not** `required='true'` and its default `value=''` is *below* its own
`min='0.50'`, so with `Advanced → Yes` the field can be submitted empty. It then reaches R as the
character `""` (`parse_args` leaves it alone: `as.numeric("")` is `NA`, so the converter declines it —
verified). And `pNbNaN <= ""` is `FALSE` for every value, **silently, with no error**:

```
pNbNaN           = 0, 0.25, 0.5, 1
LimNan = ""   -> labels: No_MV MNAR MNAR MNAR     (vs No_MV MAR MNAR with LimNan = 0.25)
c(0, 0.5, 1) <= ""  -> FALSE FALSE FALSE           (no error; the job continues)
```

So the entire dataset is classified MNAR and imputed from the global LOD, while the job log prints
"Job running in ADVANCED mode" and nothing indicates the user's threshold was ignored. `LimNan1` is
unaffected (`value='0.25'`).

**Fix:** `<param name='LimNan2' type='float' value='0.75' min='0.5' max='1' required='true' …>` (pick a
sensible in-range default rather than an empty field), and add a hard validation in the wrapper:

```r
LimNan <- suppressWarnings(as.numeric(para$LimNan))
if (length(LimNan) != 1L || is.na(LimNan) || LimNan < 0 || LimNan > 1) {
    stop("The missing value rate threshold must be a single number in [0, 1]; got: ",
         format(para$LimNan), call. = FALSE)
}
```

Add a test with `Advanced|risky|advanced=TRUE` (currently the advanced path is untested in every
combination).

---

---

### P2 — robustness, consistency, and reproducibility contract

---

#### R9 · Parameter-name drift: `mar_methods`/`mar_method`, `mnar_methods`/`mnar_method`

C19 found `para$mnar_method`; the same drift exists for MAR. Cross-referencing every `--flag` in
`<command>` against every `para$` key in the wrapper:

| name | XML flag | wrapper key used | function parameter |
|---|---|---|---|
| MAR method | `--mar_methods` | `para$mar_methods` | `mar_method` |
| MNAR method | `--mnar_methods` | `para$mnar_methods` **and** `para$mnar_method` (partial match) | `mnar_method` |

**Fix:** pick one spelling end to end. Verified there are no other orphan flags or unused keys — the
only genuine mismatch in the whole tool is `mnar_method`. `--tool_dir` is a third name in the same
table: a plumbing flag, not a tool parameter, yet it reaches `para` and so gets printed as if it were
one. R12 proposes dropping it entirely; if it is kept, at least exclude it from the printed parameter
dump.

---

#### R10 · `No_MV` features ride through the MAR imputation loop and consume the RNG stream

Features classified `No_MV` (`pNbNaN == 0`) satisfy `pNbNaN <= LimNan`, so they land in the **MAR
subset** (`:226`) and are swept through the MAR imputation loop rather than being copied through
untouched. Checked whether this can drop rows from the output: it cannot — `pNbNaN <= LimNan` and
`pNbNaN > LimNan` are exact complements, tested for `LimNan` = 0.5, 0 and `""`, and all three rows are
accounted for in every case (`nrow(MAR) + nrow(MNAR) == nrow(input)` always). So row *count* is safe.

What is *not* safe is the RNG stream: the assignment RHS is evaluated for **every** `(row, col)` pair
regardless of NA status (measured in §3.3), so complete features — and complete cells of incomplete
features — consume `mean()`/`runif()` evaluations. Consequence: **a feature's imputed values depend on
how many unrelated complete features happen to share its group subset, and on their order.** Two users
whose data differ only in features that need no imputation get different imputed values for the
features that do. That is a reproducibility trap, and it disappears once position values are computed
vectorised (R18, §3.3).

**Fix:** select rows by the computed `MARMNAR` label (C10), pass `No_MV` rows through untouched outside
the imputation loops, and assert per group:

```r
stopifnot(sum(c(No_MV = n_nmv, MAR = n_mar, MNAR = n_mnar)) == nrow(group_matrix))
```

---

#### R11 · `informations[1:12, ]` prints a wall of `NA` rows

`ImputeMSmv_function.R:203-210`: the condition is `if (nrow(informations) < 12)`, so the `else` branch
runs whenever there are **≥ 12 groups** — and it slices `informations[1:12, ]`, which is fine. But the
companion condition at `:120`, `if (length(tUclass[1, ]) < 12)`, uses `length()` of a transposed
one-row data frame — the number of groups — while indexing `[1, ]`, i.e. the *first sample's* group
label, which happens to be correct only because `t(unique(...))` has one row. Two near-identical
"12" checks with different shapes is a maintenance trap, and the `<`/`>=` boundary at exactly 12 means
a 12-group dataset prints nothing but the header. Measured for the `1:12` slice when the table is
short (defensive check): `data.frame[1:12, ]` on 5 rows silently returns 5 rows plus 7 `NA` rows — so
if the guard and the slice ever drift apart, the job log fills with `NA NA NA`.

**Fix:** one helper, `head()` instead of `[1:12, ]`, and compare against `ncol(tUclass)` rather than
`length(tUclass[1, ])`.

---

#### R12 · Hand-rolled `--tool_dir` plumbing re-implements `W4MRUtils::source_local`

*Corrected after maintainer feedback.* An earlier version of this finding claimed that
`<required_files>` "loads" the R files and therefore duplicated the `source()` call. **That was
wrong.** `<required_files>` does not source or import anything. Per the Galaxy tool-integration docs,
it "define[s] files that must be shipped from the tool directory for the tool to function properly" —
a **file-staging declaration for remote execution**, and nothing more:

> "Pulsar hacks to implicitly find referenced files from the tool directory will be disabled when this
> block is used."

Consequences, all verified:

- `ImputeMSmv.xml:26-29` lists **both** `ImputeMSmv_function.R` *and* `ImputeMSmv_wrapper.R`. That is
  **correct and should stay.** The wrapper needs to be declared because the command line invokes it as
  `Rscript '$__tool_directory__/ImputeMSmv_wrapper.R'` (`:37`) — on a Pulsar destination the script
  must be transferred before it can be run. Listing the entry point is also what the repo does
  elsewhere (`w4mconcatenate/concat.xml:8-12` lists `concatenation_wrapper.R` + 2 helpers;
  `nmr_preprocessing/NmrPreprocessing.xml:16-20` lists `NmrPreprocessing_wrapper.R` + 2 helpers).
  Because the block *disables* implicit discovery, listing every file the job touches is exactly the
  obligation the author has met — a narrower block here would be a **bug on remote destinations**, so
  do not "clean up" the wrapper entry.
- Removing the `source()` is therefore **not** an option: nothing else loads the function file. Any
  earlier suggestion to "rely on `<required_files>`" instead of sourcing must be disregarded.
- The help images (`.. image:: imputationMethod.png`, `:309`/`:377`) are served by Galaxy's static-file
  route, not staged with the job, so they need no `<required_files>` entry. Both sit under
  `static/images/` as the docs expect, matching `checkformat`, `xcms` and the rest. No action.

What *is* left is a smaller, purely maintainability issue: `ImputeMSmv_wrapper.R:26-42` hand-rolls
`get_arg()` purely to recover `--tool_dir`, then calls `base::source(file.path(tool_dir, …))`. That is
a private re-implementation of two exported `W4MRUtils` helpers the tool already depends on:

| hand-rolled | upstream equivalent (exported in the pinned 1.2.1) |
|---|---|
| `get_arg("--tool_dir")` + the `stop()` guard | `W4MRUtils::get_base_dir()` — `dirname` of `--file=`, i.e. the script's own directory |
| `source(file.path(tool_dir, "ImputeMSmv_function.R"))` | `W4MRUtils::source_local("ImputeMSmv_function.R")` |

`w4mconcatenate/concatenation_wrapper.R:9-10` does exactly that (`source_local("concatenation.R")`),
and `nmr_preprocessing/NmrPreprocessing_wrapper.R:26` does the equivalent with its own `base_dir`. So
this tool is the outlier in the repo, not the pattern.

Second-order effect of passing `--tool_dir` as an ordinary flag: `parse_args` collects every
`--flag value` pair after `--args`, so `tool_dir` lands in `para` and is dumped to the job log by
`print(para)` (`:55`) — the wrapper's own install path on the Galaxy server appears in a user-visible
log, and the flag is indistinguishable from real tool parameters in `R18`'s name-drift audit.

**Fix (small, independent of the argument-parser question):** replace the whole block at
`ImputeMSmv_wrapper.R:26-42` with

```r
source_local("ImputeMSmv_function.R")
```

(`args` is referenced nowhere else in the wrapper — `:47` re-reads `commandArgs()` itself — so it goes
too.) Delete the `--tool_dir` line from the command block (`ImputeMSmv.xml:39`). Leave
`<required_files>` untouched.

Why this is safe: `source_local` resolves relative paths against `get_base_dir()`, which is
`dirname()` of the `--file=` entry in `commandArgs(trailingOnly = FALSE)`. Measured under
`Rscript /path/to/script.R --tool_dir /opt/galaxy/tool …`, R emits
`… --file=/path/to/script.R --args --tool_dir …`, so `--file=` is present and the base dir is the tool
directory — identical to today's `$__tool_directory__`, on both local and Pulsar runners. `--args` is
inserted by R itself, which is also why `parse_args` finds the flags at all. The manual `stop()` guard
becomes unnecessary: `source()` on a missing file errors loudly by itself, with a better message than
"tool_dir not provided".

*Caveat worth knowing, not a blocker:* `get_base_dir()` falls back to `getwd()` if it finds no
`--file=` argument, so a hypothetical runner without one would silently source from the job directory
rather than fail. Under `Rscript <abs-path>` that cannot happen (measured above), and every sibling
wrapper in the repo has the same property, so accept the fallback.

Separately, and unchanged from before: `parse_args` self-deprecates —

```
message("Please, use the 'optparse' library instead of the 'parse_args' function in the future.")
```

(`W4MRUtils/R/miniTools.R:163-165`) — and 4 repo tools already use `optparse` (`batch_dispersion`,
`parsec`, `ms2snoop`, `nmr_preprocessing`). That migration is §4 **Q2**, and should not share a commit
with the correctness fixes. The `source_local` swap above is *not* part of that decision and can go in
whenever.

---

#### R13 · `rm(list = ls())` at the top of the wrapper

`ImputeMSmv_wrapper.R:18`. In a Galaxy job this is a no-op (fresh session) and it is actively hostile
when a developer sources the wrapper interactively. Drop it.

---

#### R14 · `digits` detection is wrong in both directions

`ImputeMSmv_wrapper.R:121-128` re-reads the whole data matrix as text and takes the maximum count of
characters after a `.`. Measured on the PR's own fixture:

```
digits = 12        # from raw cell "775.302795410156" — that is float artefacts, not user precision
                   # 24 002 of 70 725 cells have more than 6 decimals
=> round(x, 12) is effectively NO rounding at all
```

and for an all-integer matrix (no decimal point anywhere):

```
stripped_nb -> all NA_character_;  nchar(NA_character_) -> NA
max(..., na.rm = TRUE) -> -Inf     <- the `digits == -Inf` branch catches this, good
```

but the `-Inf` branch sets `digits <- 0`, and the *sentinel* case is subtler: `nchar(NA, allowNA=TRUE)`
is also `NA`, so the `!is.na()` filter the author probably intended is unavailable — meaning a matrix
with **some** decimal values and **some** integers reports the max over the decimals, silently
discarding precision for the integers. Also: `read.table(..., row.names = 1)` drops the identifier
column (fine) but the whole matrix is read twice (once by `import3`, once here) purely to derive a
rounding parameter.

**Fix:** this is a hidden, undocumented behaviour with a user-visible effect on output precision.
Either (a) expose it as an explicit tool parameter with a documented default, or (b) drop rounding
entirely and let R's default 15-significant-digit output stand, or (c) if keeping detection, base it
on a precision *cap* rather than observed digit counts and log why. Given "rounding precision applied
during imputation" (`:63`) is not mentioned in the user-facing help at all, (a) or (b) — decide with
the maintainer (§4 **Q1**). Whatever is chosen, note that it changes the output fixtures (§3.1).

---

#### R15 · A numeric-looking grouping column name silently selects the wrong column

`W4MRUtils::parse_args` (`miniTools.R:156-200`) applies a boolean converter then a numeric converter
to **every** value, first success wins (`convert_parameters`, `miniTools.R:230-246`). So `--col 1`
for a sample-metadata column literally named `1` arrives as `para$col = 1` (numeric). Verified:

```r
p$col <- 1;  p$col %in% colnames(s)          # TRUE   (R coerces for the comparison)
colnames(s[1])                                # "1"    <- right only because it IS the first column
# with the group column NOT first:  s[[1]] name == "z"  instead of "1"   <- silently wrong column
```

The guard at `ImputeMSmv_wrapper.R:141-145` passes, then `s[[col]]`/`get(col)` resolve by position.
**Fix:** in the validator, `stopifnot(is.character(col))` and coerce explicitly:
`col <- as.character(para$col)`; then look the name up. Note `advanced`/`setSeed` are **not** affected
— `convert_parameters` tries the boolean converter first and stops on success, so `"TRUE"` becomes
`TRUE` (verified), which is why `para$advanced == TRUE` at `:150` works. Prefer `isTRUE()` anyway.

---

#### R16 · `assign()`/`get()`/`mget()` with group-derived variable names

`ImputeMSmv_function.R:132-141, 146-199, 220-248, 253-311, 318-332` build object names from user data
(`paste("DataMatrix_", tUclass[i], sep="")`) and read them back with `get`/`mget`. Group labels are
free text, so **suffix collisions silently corrupt the output**. Reproduced:

```
groups "a" and "a_MV"   (both legal W4M labels; the help only discourages - + * / é)
  group "a"   creates DataMatrix_a  and  DataMatrix_a_MV   (its MV table)
  group "a_MV" creates DataMatrix_a_MV                     <- COLLIDES

no error raised; output data frame gains duplicate column names and the bookkeeping columns leak:
      a1   a2   b1   b2 NA   NA        <- duplicated NA col names, and 1/0.25 bookkeeping values
f1  1.00 1.09 1.00 0.75  1  0.25
dup colnames: TRUE
```

The XML help (`:413`) says "It is recommended not to use special characters" — insufficient, because
`_` is explicitly *allowed* and `_MV`/`_MAR`/`_MNAR`/`_noMAR`/`_noMV`/`_MARMNAR` are exactly the
suffixes the code appends, so a label like `Cond_MV` or `QC_noMV` is a self-inflicted corruption.

**Fix:** remove the pattern entirely — keep per-group results in a `list()` keyed by group name (see
§3). No name sanitisation is needed once no identifier is synthesised from data. If a transitional
guard is wanted, `check_name <- make.names(...)` plus a `setequal()` assertion is not enough; the
list-based rewrite is the fix.

---

#### R17 · `tUclass` derived via `t(unique(sampleMetadata[colClass]))`

`ImputeMSmv_function.R:114`. Works (verified: the transposed `data.frame` becomes a character matrix
and group labels compare correctly, including when the column is a factor — `t()` forces
`as.character`, so no factor-code confusion), but it is obscure, it silently loses duplicate-free
guarantees if the subset ever returns >1 row, and it forces the awkward `tUclass[i]`/`tUclass[1, ]`
indexing seen in R11. **Fix:** `groups <- unique(as.character(sampleMetadata[[colClass]]))`, used as a
plain character vector, plus `length(groups) > 0` validation.

---

#### R18 · MAR/MNAR position values are recomputed from the matrix being mutated in place

The imputation loops read their reference value from the same object they write into
(`:234`/`:240` `mean(Mat_noMar_temp_2[seq(2, nrow(...)), j], na.rm = TRUE)`), so each imputation
consumes the values imputed before it. Measured, default `SDmean = 0.1`:

```
feature with observed 100 (×4) and 4 MVs, expected band 100 * U[0.9, 1.1] = [90, 110]
  imputed: 99.000  108.782  100.284  111.268      <- 111.27 outside the documented noise band
with SDmean = 0 (noise exactly 1): 100 100 100 100 (the feedback loop is masked only when noise == 1)
median variant, observed 100,104 -> median 102, expected 102 ±10% = [91.8, 112.2]
  imputed: 98.94  93.00                            <- second draw pulled by the first
```

The mean case is self-stabilising (inserting the mean does not move the mean), which is why the noise
happens to be the only visible symptom; the **median** case is genuinely order-dependent, and the
documented behaviour ("the mean or median of the group replicates", `ImputeMSmv.xml:108`, `:279-282`)
does not describe position values that drift.

**Fix:** compute the per-(feature, group) position value **once**, from the observed values only,
before the imputation sweep; then apply noise per cell. This removes the feedback loop, and it is
also what makes the vectorised rewrite in §3 straightforward.

---

#### R19 · `pmin = 0` and `percLOD` extremes are accepted and produce degenerate data

`ImputeMSmv.xml:131` allows `min='0'` for `percMIN`. With `pmin = 0` the tool imputes exact zeros
(verified: matrix becomes `100 0 0 / 0 200 0`), which then flow into `plt_density`'s `log10` —
`log10(0) = -Inf`. `density()` tolerates it here, but downstream W4M tools that log-transform will
break, and imputing literal 0 is precisely what the tool's own help says it exists to avoid
(`:373`, "It replaces the common practice of setting missing values to zero"). **Proposed fix:** set
`min` to a small positive value (e.g. `0.01`) for `percMIN`, and validate `percLOD ∈ (0, 1]`.
`percLOD = 0` yields the matrix minimum as LOD — arguably intended, so document rather than reject.

---

#### R20 · Everything-NA input yields "successfully imputed" `NA`s

Reproduced: an all-`NA` data matrix (the wrapper's `sum(is.na(d)) == 0` check at `:102` passes, since
there *are* missing values) runs to completion, prints the LOD banner, and writes a matrix that is
still entirely `NA` (`"f1"\tNA\tNA`). Galaxy marks it OK. **Fix:** the validator (C8) should reject
features with no observed value anywhere (R4/C17), and the wrapper should additionally verify
`sum(is.na(result$Datamatrix)) == 0` before writing — a cheap post-condition that turns any residual
silent failure into a job error.

#### R21 · The RNG stream is consumed by every cell, so imputed values depend on loop order and on unrelated features

Verified: the assignment `x[k,j][is.na(x[k,j])] <- round(ref * round(runif(1, …), 2), digits)`
evaluates its **RHS unconditionally** — R evaluates the right-hand side of an assignment before
deciding what to overwrite. A matrix with **no** NA at all still triggers `nrow × ncol` evaluations:

```
2x2 matrix, 0 NAs -> 4 runif-equivalent evaluations
```

Two consequences:

1. **The two loops walk the matrix in opposite orders.** The MAR branch transposes first (`:227`) and
   then iterates `k` over *samples* and `j` over *features* — i.e. samples outer, features inner — while
   the MNAR branch (`:283`, `:301`) iterates features outer, samples inner. So the same `(feature, cell)`
   pair receives a different noise draw depending on which branch it falls into, and on the order of
   features in the file. Measured, flipping the iteration order on a 4×6 block changes **all 10** imputed
   values (`identical? FALSE`), even though no value is written to a non-NA cell.
2. **A feature's imputed values depend on unrelated features.** Because every cell consumes a draw, the
   position of feature *f*'s NAs in the stream is determined by how many features precede it and how
   many *complete* features are in the subset — so deleting a feature that needs no imputation changes
   the imputed values of the features that do. This makes the recorded-seed reproducibility promise
   fragile in a way users cannot predict.

Not a numeric-correctness bug (the values drawn are valid), but it is a reproducibility defect and it
silently pins the implementation to its exact iteration order — any "harmless" refactor, including the
one proposed in §3.3, changes the output. That must be a stated, deliberate change rather than a
surprise discovered when the fixtures stop matching.

**Fix:** compute each feature's position value once (R18), then draw noise only for the cells actually
being imputed — one draw per NA, in an order determined by the data matrix's own row/column order
(`impute_na()` in §3.3). Then assert, in a test, that the result is invariant to the order in which
*unrelated* features appear:

```r
# same NA pattern, features reordered -> each feature must impute to the same value
stopifnot(isTRUE(all.equal(res$Datamatrix[f, ], res_perm$Datamatrix[f, ])))
```

---

### P3 — style, lint, documentation

---

#### R22 · roxygen `@return` block is self-contradictory

`ImputeMSmv_function.R:65-69` has one `@returns` ("A list containing 3 W4M-compatible tables") and
**three** `@return` tags describing the three elements as if they were separate return values.
roxygen2 keeps only one. Fix: one `@return A named list with elements \code{Datamatrix},
\code{sampleMetadata}, \code{VariableMetadata}`, and document each element under `@param`/`@details`.

---

#### R23 · Typographic and spelling defects that will reach users

- `ImputeMSmv_function.R:1`-style comment headers, and the XML help, use en/em dashes and curly
  quotes in identifiers-adjacent and `cat()` contexts: `spectrometry–based` (`:16`, `ImputeMSmv.xml:231`),
  `"` `"` — these render as mojibake in some Galaxy log viewers and break `grep`-ability. `.shed.yml:4`
  also contains `–`.
- Typos in user-facing help: `recommanded` (`:92`), `indenpendently` (`:253`), `remcommended` (`:262`),
  `defaut` (`:293`, `:447`), `strutured` (`:310`, `:378`), `propotion` (`:149`, `:455`), `positionned`
  (`:371`), `additionam` (`ImputeMSmv_function.R:33`), `valu statitics` (`:116`), `reproductibility`
  (`:96`, `:179`), `Comapre` (`:377`), `matix` (`:439`), `Chosen` (`:450`), `Na ->` (`:451`). Also
  `Paramters set:` (`ImputeMSmv_wrapper.R:54`).
- `ImputeMSmv_wrapper.R:51` greets the user with **"Welcome to MissMDA, a Galaxy tool dedicated to the
  handling of missing values from Mass Spectrometry data."** — the first line of every job log. The name
  is wrong, and it is not merely a misspelling of this tool: **missMDA is a real, unrelated package**
  (Josse & Husson, *J. Stat. Soft.* 70(1), 2016, DOI 10.18637/jss.v070.i01), which imputes under an
  MCAR assumption using PCA-family/iterative-SVD methods — a different method family from this tool's
  MAR/MNAR classification, and unrelated to the method paper actually cited here (Delaporte, Cladière &
  Camel 2019, *Chemom. Intell. Lab. Syst.* 188:54‑62 — see C2). Left as is, it tells users they are
  running a different method than they are, and anyone searching the instance for "missMDA" lands on the
  wrong tool. Should read ImputeMSmv, and the sentence should describe the MAR/MNAR approach.
  *Checked:* a repo-wide grep finds `missMDA` nowhere else in this branch, so this is a one-line fix
  local to this wrapper, not a pattern copied from a sibling tool.
- Function header comment says `Tool: ImputeMsmv` (`ImputeMSmv_function.R:2`) — letters transposed vs
  `ImputeMSmv`.

**Fix:** run a spell pass over the help and the `cat()` strings; replace non-ASCII punctuation with
ASCII in code and log output (keep accents in the human credits block, where they are correct); rewrite
the welcome line. If the MissMDA lineage is intentional as prior art, cite it properly in
`<citations>` rather than in the greeting.

---

#### R24 · Dead / inconsistent code

- `:243-245` `else stop("Non valid method…")` for `mar_method` is unreachable — `match.arg` already
  rejected anything else at `:93`. Same for `:312-314` (`mnar_method`) and `:266-268`
  (`mnar_computation`). Keep exactly one validation point (`match.arg`); the duplicated messages will
  drift.
- `:282` and `:300` guard `if (nrow(get(...)) > 0)` but `:279`/`:298` already fetched the same object
  via `get()`; and in the `min` branch (`:301-306`) the loop runs unguarded when `nrow == 0` would be
  false — fine, but the asymmetry with the `lod` branch (`:283`) is unexplained.
- `:213` `mget(nom_DataMatrix_gpMARMNAR[seq(nom_DataMatrix_gp)])` — `seq(nom_DataMatrix_gp)` on a
  character vector is `seq_along(...)` written obscurely; same at `:332` with `seq_along`. One form,
  used consistently.
- `:228` removes bookkeeping columns by *row* name after `t()`, `:289`/`:308` remove the same columns
  by *column* name — two idioms for one operation.
- `:321-325` `rbind` → `order()` → `[, -1]` restores input order by index; but `order()` on a
  possibly-duplicated `numero` (possible right after C10's fix if classes overlap) silently reorders —
  make the order restoration explicit (`[rownames(ordre), ]`).
- `:446` `plt_matrix(df_original)`: the roxygen `@param df` (`:437`) does not match the actual argument
  name `df_original` — R CMD check / `lintr` objection, and wrong in generated docs. Same class:
  `plt_group`'s `@param ID` vs the `ID` argument is fine, but `plt_type`/`plt_pie` never use `ID` for
  anything except group extraction — verified used. Also `plt_density(df1, df2)` documents "df1 Data
  matrix before imputation" but the caller passes the *quoted* `DM_out` (`ImputeMSmv_wrapper.R:237`),
  hence the `fix_head()` call inside the plot function (`:400`) — a plot function repairing its
  arguments is a layering smell; see §3.

---

#### R25 · `plt_density` re-derives row names from the output it was given

`:400` `df2_fixed <- fix_head(df2)` exists only because the wrapper passes `DM_out` (identifier column
still attached, `ImputeMSmv_wrapper.R:237`) instead of `result_tables$Datamatrix`. Consequence: the
"after" distribution silently depends on whether column 1 happens to be an ID column, and `fix_head`
`stop()`s if the first *data* column has duplicates — a plotting-only failure mode that kills the job
after imputation already succeeded. **Fix:** pass the tidy matrix; delete `fix_head` from the plot
function.

---

#### R26 · `plt_group` divides by `n_features` without a zero guard; percentages are per feature×group cell

`:515` `col_vals <- c(No_MV, MAR, MNAR) / n_features * 100` with `n_features <- nrow(groups)`; a
zero-feature matrix gives `NaN`. Also, unlike `plt_type`, `plt_group` normalises per group by feature
count, so the two plots look comparable but use different denominators — worth one help sentence.
**Fix:** reuse the `pc()` helper from C15.

---

#### R27 · `requirements`: pin consistency

`ImputeMSmv.xml:21-24` pins `r-base 4.5.2` and `r-w4mrutils 1.2.1`. In this repo the common r-base pin
is `4.3.3` (3 tools) with `4.5.2` used by 2, and `r-w4mrutils` by `1.0.0` (2 tools) or `1.2.1` (1). The
pins are internally consistent with `README.md:33-34` ("4.5.2 recommended", "1.2.1"), so **keep 4.5.2**
— but confirm the tool actually runs on the Galaxy-pinned tool dependency (I could only execute on R
4.3.3, where all findings reproduce), and consider `@TOOL_VERSION@`-style templating if the repo moves
that way (1 tool already does).

---

## 3. Refactoring: readability, maintainability, duplication

The core observation: **`imputemsmv()` is ~270 lines that manipulate data by constructing R variable
names from user input.** Almost every finding above (R1, R2, R5, R10, R16, R18) is a downstream
consequence of that one design choice, and none of the fixes require new capability. The refactor is
therefore the *fix*, not a nicety.

### 3.1 Ground truth to preserve

Verified before touching anything: the committed fixtures are **exactly reproducible** from the
current code. Running the pipeline as the wrapper does, with `lim_nan=0.5`, `seed=123`,
`mnar_method="lod"`, `percLOD=0.001`, `lod_computation="mean"`, `mar_method="mean"`,
`sd_mar=0.1`, `sd_mnar=0.2`, `advanced=FALSE`, `digits=<detected 12>`:

```
DM fixture: IDENTICAL
SM fixture: IDENTICAL
VM fixture: IDENTICAL
```

So this is the regression oracle. **Any refactor that does not reproduce these three files byte for
byte has changed behaviour and must be justified explicitly.** (This will be *deliberately* broken by
C20 `quote=FALSE`, by R14 if rounding changes, and by C10/C12 once the classification/imputation
actually follows the documented method — each of those needs a fixture regeneration with a note in the
PR describing the scientific change.)

### 3.2 Replace `assign`/`get`/`mget` with lists (removes R16; enables R1, R2, R10)

Current: 5 parallel `nom_*` name vectors, 6 `assign()` sites, 11 `get()`/`mget()` sites, and names
synthesised from group labels. Replacement:

```r
split_by_group <- function(DataMatrix, sampleMetadata, groups) {
    lapply(groups, function(g) {
        samples <- rownames(sampleMetadata)[as.character(sampleMetadata[[colClass]]) == g]
        DataMatrix[, samples, drop = FALSE]
    }) |> setNames(groups)
}
```

then `lapply()` over the list for classification, MAR imputation, MNAR imputation, and reassembly with
`do.call(rbind, ...)` / explicit row reindexing. This deletes `ordre`/`numero`, the `[, -1]` index
column dance, and the entire class of name collisions. It also makes the MAR and MNAR steps unit-testable
in isolation, which they currently are not (see §3.5).

### 3.3 Remove the three near-identical imputation loops (removes R18, R5; the largest duplication)

`:232-236` and `:238-242` differ in exactly one token (`mean` vs `median`). `:283-287` and
`:301-306` differ in how the reference is computed. All four are the same nested `for (row) for (col)`
with `x[k,j][is.na(x[k,j])] <- round(ref * round(runif(1, 1-sd, 1+sd), 2), digits)` (verified that the
double-assignment idiom works, but it is needlessly clever — `x[k,j]` is already a scalar).

One vectorised primitive replaces all of them:

```r
#' Replace NA in `x` by `ref * noise`, where noise ~ U(1-sd, 1+sd), rounded to `digits`.
impute_na <- function(x, ref, sd, digits) {
    na <- is.na(x)
    if (!any(na)) return(x)
    noise <- round(runif(sum(na), min = 1 - sd, max = 1 + sd), digits = 2)
    x[na] <- round(ref[na] * noise, digits = digits)
    x
}
```

with `ref` precomputed as a full matrix by broadcasting a per-feature or per-(feature,group) position
value — which is exactly the "compute the position value once" fix from R18, and removes the
`mean`/`median` branch duplication (one `FUN` argument). Note the `round(..., digits = 2)` on the
noise multiplier deserves a comment plus a test: measured, it **collapses the noise to a constant 1.00
whenever `sd < 0.005`** (sd=0.01 → 3 distinct multipliers; sd=0.005 or 0.001 → exactly one), silently
defeating the tool's stated purpose of "avoid artificial ties" (`ImputeMSmv.xml:281`). Either drop the
inner rounding or raise the parameters' `min`.

**This rewrite is not behaviour-preserving, and the plan must not pretend otherwise.** The current
loops evaluate their RHS once per `(row, col)` pair — verified: a matrix with **zero** NAs still
triggers 4 evaluations on a 2×2 matrix — so they consume `n_features × n_replicates` draws per group,
while `impute_na()` consumes one per NA. Same seed, different values:

```
matrix with 3 NAs of 6, seed 1
  nested loops (4 draws) : 1080, 970, 1010
  vectorised   (3 draws) :  950, 970, 1010      identical? FALSE
```

Consequences, both good but both needing to be said out loud:

- The **output values change** for the same seed, so the `lod`+`mean` fixtures **must be regenerated**
  in this step. §3.1's byte-identity gate therefore applies to steps 2, 5, 6 only.
- The dependency on unrelated cells disappears: after this change a feature's imputed values depend
  only on that feature's own data, which removes R10's reproducibility trap and makes T4/T5 meaningful.
- The `seed` a user recorded before the upgrade will not reproduce the old run. Worth one line in the
  tool help ("imputation results for a given seed changed in version 1.0.1"), since reproducibility is
  the explicit promise of the `Reproducibility` section (`ImputeMSmv.xml:96`).

### 3.4 Single source of truth for MV classes, bookkeeping columns, and percentages

- One `classify_mv(pNbNaN, rep_num, LimNan, advanced)` returning the label vector, used for the
  variable-metadata output **and** for row selection (fixes C10 by construction, R10).
- One `BOOKKEEPING <- c("NbNaN", "pNbNaN")` constant and one `drop_bookkeeping()` helper (replaces
  `:228`, `:289`, `:308`).
- One `pc()` percentage helper + one fixed `MV_COLORS` map, used by `plt_group`, `plt_type`, `plt_pie`
  (fixes C15, C16, R26; the three plot functions currently each re-derive counts from the same
  variable metadata with three different idioms — `for`+`cbind`, `colSums`, `sum` — and two different
  colour sources, hard-coded triples and `hcl.colors`).
- Validation in one wrapper-side function (C8), not six blocks in two files.

### 3.5 Test coverage (currently 3 test blocks for 2×3×3×2×5 parameter space)

`<tests>` has 2 real cases (1 pass, 1 `expect_failure`), all on `lod`+`mean`; the PDF output
`plot_out` is **never asserted** by any test. Untested and measurably broken: `min` (R3),
`lod_computation=median`, `mar_method=median`, `advanced=TRUE`, single-MAR-feature groups (R5),
interleaved groups (R1), VM/SM in a different order (R2), unblocked designs, the no-replicate design
the help recommends (R5), and the all-NA guard (C17). Proposed additions, each independently valuable:

| # | Test | Asserts |
|---|---|---|
| T1 | `mnar_method=min` | imputed values within `[pmin*min(obs), pmin*max(obs)*(1+SDLOD)]`; no `Inf`/`NaN` |
| T2 | `advanced=TRUE`, 2-replicate group | MV class in VM == strategy actually used (regression for C10) |
| T3 | group with exactly 1 MAR feature | no crash (regression for R5) |
| T4 | interleaved groups, 1 feature per class | every observed input cell unchanged in the output column of the same name (regression for R1) |
| T5 | VM/SM rows in reversed order | `mv_type_*`/`MV_percentage` attached to the right ids (regression for R2) |
| T6 | one feature 100 % NA | `expect_failure='true'` with the feature id in the message (regression for C17/R4) |
| T7 | no-replicate design (one sample per group) | no crash, or documented failure (regression for R5) |
| T8 | `plots` = all selected / none selected | PDF output exists and is non-empty (regression for C5) |
| T9 | `mar_method=median`, `lod_computation=median` | golden outputs (these paths never execute in CI today) |

T4 and T5 are the two that would have caught the P0s. `compare="size"` or a checksum is acceptable for
the PDF; planemo also supports `change_exe`/`assert_content` for content checks.

### 3.6 Suggested commit sequence

Small, independently testable, correctness first. Each step's validation is the §3.1 fixture
comparison plus its own new test.

| # | Commit | Files | Risks | Validation |
|---|---|---|---|---|
| 1 | Move all input validation into the wrapper validator; fix the inert all-NA guard | `ImputeMSmv_wrapper.R`, `.xml` (2 new `expect_failure` tests) | Rejects inputs that previously "succeeded" (R20) — intentional; must be in the release note | T6; fixtures unchanged |
| 2 | Reindex outputs by name (DM columns, VM rows, SM rows) instead of positional rename | `ImputeMSmv_function.R` | None on blocked data; fixes interleaved designs | T4, T5; fixtures byte-identical |
| 3 | Single MV classification source; select MAR/MNAR rows by label | `ImputeMSmv_function.R` | **Changes imputed values** on small-group data (scientific change) | T2; regenerate DM/VM fixtures with a written justification |
| 4 | Fix the `min` reference (per-feature, all samples, no bookkeeping cols) | `ImputeMSmv_function.R` | **Changes imputed values** for `min`; `Inf` risk if step 1 not merged first | T1; new `min` fixture |
| 5 | Fix `seq_len(NULL)` collapse; `drop = FALSE` throughout | `ImputeMSmv_function.R` | Low | T3, T7 |
| 6 | Reproducibility: seed in exactly one place; clock seed otherwise; reconcile XML help wording | `.xml`, `_wrapper.R`, `_function.R` | Low | rerun T1–T5 twice, identical |
| 7 | Remove function-argument defaults (C6) + resolve `mnar_computation`/`percLOD`/`pmin` plumbing | `_function.R`, `.xml` | **Must follow step 6's command-block change**, else `min` breaks (R6) | T1 with `min` |
| 8 | Refactor: lists instead of `assign`/`get`; vectorised `impute_na`; position values computed once | `_function.R` | **Changes imputed values for a given seed** (RNG consumption order, see §3.3 and R21) — cannot keep fixtures byte-identical | regenerate all 3 fixtures; assert every non-NA input cell unchanged; assert imputed values inside the documented noise band |
| 9 | Plot fixes: `pc()`/`MV_COLORS`, pie mask, empty-PDF page, `plt_density` argument, `on.exit` inside a function | `_function.R`, `_wrapper.R`, `.xml` | Low | T8 |
| 10 | Named result access (C4), `quote=FALSE` for DM (C20), `detect_errors` (C3), `para$mnar_methods` (C19), drop `rm(list=ls())` (R13), `source_local` swap + drop `--tool_dir` (R12) | `wrapper.R`, `.xml` | fixture header changes; R12 needs one real tool test run, not just unit tests | regenerate fixtures; tool tests pass |
| 11 | Metadata/docs: `.shed.yml` (C1, C2), README (C18 + provenance note), `<creator>` + drop Contact (C9), typos/ASCII (R23), roxygen (R22), `digits` decision (R14) | `.shed.yml`, `README.md`, `.xml` | none | `planemo shed_lint`, `planemo lint --xml_schema` |
| 12 | Optional / follow-up: `optparse` migration (R12's remaining half, gated on Q2), `param` renames, test matrix expansion (T9) | all | touches every flag | full tool-test run |

Steps 1–5 are the ones that stop wrong science from being produced. Steps 8 and 12 can be follow-up
PRs if the contributor prefers, but **8 should not be deferred past a second release** — the `assign`
pattern is the root cause of four findings.

---

## 4. Questions needing the PR author / maintainers

These change the fix, so they are not mine to assume:

- **Q1 (R14, `digits`).** Is auto-detecting output precision from the input file's decimal notation
  intended behaviour? It currently yields 12 for the test fixture (float artefacts → no rounding).
  Options: explicit parameter / no rounding / keep with docs. Affects the fixtures.
- **Q2 (R12, optparse half only).** Migrate to `optparse` now (repo direction, `parse_args`
  self-deprecates) or keep `parse_args` for this PR and do it separately? The `source_local` half of
  R12 is *not* gated on this answer and can be done either way.
- **Q3 (C10 + R18).** Confirm the intended semantics for a group with `< 3` replicates: is
  "any missingness ⇒ MNAR, imputed with the MNAR strategy" correct, i.e. should the small-group rule
  override `LimNan` for **both** label and strategy? My assumption throughout is yes.
- **Q4 (R19).** Should `percMIN = 0` be permitted? The tool's stated purpose is to avoid zero-filling.
- **Q5 (R27).** Confirm `r-base 4.5.2` + `r-w4mrutils 1.2.1` are available as Galaxy dependencies, and
  that the tool is expected to work on the repo's more common `4.3.3`.
- **Q6 (C9).** Which email/support address should `<creator>` carry? Do not want to invent one.
- **Q7 (C18, raised while implementing step 11).** Which dataset are `inDmShT1` / `inSmShT1` /
  `inVmShT1` derived from? `README.md` says "derived from publicly archived datasets available on
  SoftwareHeritage" with no accession, and the Software Heritage link in the help archives the
  original R *code*, not data, so it cannot be the source. This matters more here than usually: the
  tool test compares all three outputs byte-for-byte, so an unrecorded fixture source is an
  unreproducible test. (These three names appear in no other tool in this repo, so they are not
  borrowed from a sibling's test-data.)
- **Q8 (C18, withdrawn after measuring).** I first wrote that this README sentence was simply wrong.
  It is not. `test-data` is 2.75 MiB in total, but *every individual file* is under 1 MB (largest:
  `outDmShT1` at 979 018 bytes), which is the natural reading of "keeps the test dataset files under
  1MB". The inputs alone total 1.29 MiB, so the sentence is only true file-by-file — at most worth an
  "each", not a correction.

---

## 5. Verification checklist before the PR is ready

```bash
cd tools/imputemsmv
planemo lint .                                   # + --xml_schema
planemo test  .                                  # all blocks, incl. the new T1–T8
planemo shed_lint --tools --ensure_hg1_consistency .
Rscript -e 'lintr::lint_dir(".")'                # CI job "Lint R scripts"
```

- [ ] `planemo test` green with T1–T8 added (currently 3 blocks; target ≥ 10)
- [ ] fixtures byte-identical after steps 2, 5, 6 (behaviour-preserving)
- [ ] fixtures regenerated **and the seed-semantics change announced** after step 8 (RNG order, R21)
- [ ] fixtures regenerated **with a justification paragraph** after steps 3, 4, 10 (behaviour-changing)
- [ ] every reviewer comment either fixed or answered with the measured behaviour
      (in particular C11 → "refuted, `quantile()` flattens data frames"; C13 → "top-level `on.exit`
      is a no-op on R 4.3.3, fixed as hygiene, it did not abort the job")
- [ ] `README.md` and `.shed.yml` agree with the XML on id, version, description
- [ ] all four outputs asserted by at least one test, including `plot_out`
- [ ] XML help updated wherever behaviour changed (small-group rule, seed semantics, `percMIN` range,
      plots-when-none-selected) — and the "no replicates" guidance either fixed or made to work (R5)
- [ ] non-ASCII punctuation removed from code and `cat()` output (R23)

---

## Appendix A — reviewer-comment disposition at a glance

| ID | Author | Location | Verdict | Severity | Fix |
|---|---|---|---|---|---|
| C1 | hechth | `.shed.yml:2` | accepted | P3 | §1.1 |
| C2 | hechth | `.shed.yml:5` | accepted (DOI is the method paper — note it) | P3 | §1.1 |
| C3 | hechth | `xml:33` | accepted | P3 | §1.1 |
| C4 | hechth | `wrapper:207` | accepted | P3 | §1.1 |
| C5 | hechth | `xml:165` | accepted; current answer = 0-page PDF | P2 | §1.1 |
| C6 | hechth | `function:88` | accepted; required for the `min` path | P2 | §1.1, R6 |
| C7 | hechth | `function:101` | accepted (conclusion; mechanism refined) | P2 | §1.1, C14 |
| C8 | hechth | `function:110` | accepted | P2 | §1.1 |
| C9 | lecorguille | `xml:219` | accepted | P3 | §1.2 |
| C10 | Copilot | `function:226` | **re-scoped**: real bug, wrong mechanism/example | **P1** | §1.3 |
| C11 | Copilot | `function:259` | **REFUTED** | — | §1.3, reply only |
| C12 | Copilot | `function:302` | confirmed (both halves) | **P0** | §1.3, R3 |
| C13 | Copilot | `wrapper:221` | mechanism refuted; fix anyway | P3 | §1.3 |
| C14 | Copilot | `function:102` | confirmed (= C7) | P2 | §1.3 |
| C15 | Copilot | `function:590` | confirmed; draws literal `NaN%` | P2 | §1.3 |
| C16 | Copilot | `function:661` | confirmed; symptom differs (colours, not labels) | P2 | §1.3 |
| C17 | Copilot | `wrapper:115` | confirmed; guard is inert on all inputs, fixtures hit it | **P0/P1** | §1.3, R4 |
| C18 | Copilot | `README:8` | confirmed | P3 | §1.3 |
| C19 | Copilot (overview) | `wrapper:164` | confirmed; works via partial `$` matching | P2 | §1.3, R9 |
| C20 | Copilot (overview) | `wrapper:210` | confirmed; inconsistent with SM/VM and repo norm | P3 | §1.3 |

## Appendix B — new findings index

| ID | Severity | One-line summary (taken from the finding's own heading) | Location |
|---|---|---|---|
| R1 | **P0** | Output columns are mislabelled whenever replicate groups are not contiguous | `function:332-333` |
| R2 | **P0** | Variable metadata and sample metadata are `cbind`ed positionally | `function:336, 339` |
| R3 | **P0** | The `min` MNAR method is broken (C12) — see C12 for the fix | `function:302` |
| R4 | **P0** | The all-NA-feature guard is inert (C17) — see C17 for the fix | `wrapper:114-118` |
| R5 | **P1** | Crash when any group has exactly one MAR feature | `function:227-228, 232, 238` |
| R6 | **P1** | `mnar_computation` never reaches the function on the `min` path | `xml:63-68`, `wrapper:192` |
| R7 | **P1** | `percLOD`/`pmin` share the same problem | `xml:63-68` |
| R8 | **P1** | `LimNan2` has `value=''` and no `required`, silently collapsing the whole matrix to MNAR | `xml:154` |
| R9 | P2 | Parameter-name drift: `mar_methods`/`mar_method`, `mnar_methods`/`mnar_method` | `wrapper:164, 230` |
| R10 | P2 | `No_MV` features ride through the MAR imputation loop and consume the RNG stream | `function:226, 278, 297` |
| R11 | P2 | `informations[1:12, ]` prints a wall of `NA` rows | `function:120, 203-210` |
| R12 | P2 | Hand-rolled `--tool_dir` plumbing re-implements `W4MRUtils::source_local` | `wrapper:26-42`, `xml:39` |
| R13 | P2 | `rm(list = ls())` at the top of the wrapper | `wrapper:18` |
| R14 | P2 | `digits` detection is wrong in both directions | `wrapper:121-128` |
| R15 | P2 | A numeric-looking grouping column name silently selects the wrong column | `wrapper:141` |
| R16 | P2 | `assign()`/`get()`/`mget()` with group-derived variable names | `function:132-332` |
| R17 | P2 | `tUclass` derived via `t(unique(sampleMetadata[colClass]))` | `function:114` |
| R18 | P2 | MAR/MNAR position values are recomputed from the matrix being mutated in place | `function:234, 240` |
| R19 | P2 | `pmin = 0` and `percLOD` extremes are accepted and produce degenerate data | `xml:131` |
| R20 | P2 | Everything-NA input yields "successfully imputed" `NA`s | `wrapper:102, 210` |
| R21 | P2 | The RNG stream is consumed by every cell, so imputed values depend on loop order and on unrelated features | `function:227, 234, 283, 301` |
| R22 | P3 | roxygen `@return` block is self-contradictory | `function:65-69, 437, 444` |
| R23 | P3 | Typographic and spelling defects that will reach users | `wrapper:51`, `function:2`, `xml` help |
| R24 | P3 | Dead / inconsistent code | `function:243, 312, 213, 228, 289` |
| R25 | P3 | `plt_density` re-derives row names from the output it was given | `function:400` |
| R26 | P3 | `plt_group` divides by `n_features` without a zero guard; percentages are per feature×group cell | `function:505, 515` |
| R27 | P3 | `requirements`: pin consistency | `xml:21-24` |

---

## Appendix C — implementation record

What was actually implemented, and where measurement contradicted this report. Steps follow §3.6.

### Step 9 (plots) — implemented

Covered: C5, C13, C15, C16, R25, R26, plus §3.4's shared `MV_COLORS` / `pc()`.

**Four corrections to this report, all from measuring rather than reasoning:**

1. **C13's `on.exit` benefit is smaller than §1.3 implies.** I wrote that an error inside a plot
   "leaves Galaxy with a truncated PDF". Tested by patching `graphics::image` to `stop()` mid-run
   (`trace(..., where = asNamespace("graphics"))`, injected via `R_PROFILE_USER`) and running the tool
   before and after the change: **both** exit 1 and leave a *well-formed* PDF holding the 4 pages
   completed so far, because grDevices finalises open devices at process exit. Moving the plotting
   into `make_plots()` does make the handler fire on the error path — verified separately, a
   function-scoped `on.exit` runs during unwinding where the top-level one does not — but it fixes
   dead code, not an observable defect. Severity stays P3.
2. **C5 was fixed without the `<option value='None'>` this report proposed.** `None` appears nowhere
   in the XML or README — it was purely dead code in the wrapper — and a `None` checkbox alongside
   five others invites the contradiction "`None` + `DP`". Instead, an *empty* selection now emits a
   "NO PLOT SELECTED" page, the help says so, and the wrapper filters `None`/`""` out of the
   selection so a stray `None` means "no plot" and `None,DP` still plots DP. Verified end to end:
   the pre-fix code on `--plots ''` produced `PDF document, version 1.4, 0 page(s)`; the new code
   produces 1 page.
3. **C15's `NaN%` lives in `plt_type`, not `plt_group`.** This report cites `function:590`, which is
   inside `plt_group` — but `plt_group` normalises *per group* by the feature count, so its `NaN`
   needs `nrow(groups) == 0` and is unreachable in practice. The per-category `x / sum(x)` with no
   zero guard — the shape Copilot described — is in `plt_type`, and it reproduces on ordinary data
   (`MAR` absent overall ⇒ the `MAR` column is all `NaN` and is drawn). Both are now guarded:
   `plt_type` drops empty categories outright, `plt_group` is routed through `pc()`.
4. **R25's "silently dropped column" does not happen.** §2 R25 says the `fix_head()` inside
   `plt_density` exists because the wrapper passes `DM_out`, implying an identifier column was being
   mishandled. Measured against the pre-fix code: `fix_head()` correctly stripped the id column and
   the "after" curve saw all 6 intensities — identical to the tidy-matrix result. So the defect was
   layering plus a *latent* failure mode (`fix_head()` `stop()`s on duplicate or NA identifiers,
   which would kill the job after imputation had succeeded), not a wrong curve. The change is still
   worth making — the contract is now "intensities only", and an identifier column fails loudly at
   the call site instead of being rewritten inside a plotting function — but it must not be
   described as a numerical fix.

**Two deliberate visual changes, both beyond the letter of the findings:**

- `plt_pie` percentages now divide by the classified cell count instead of `features × groups`.
  Previously, when the variable metadata held any cell that was neither `No_MV`, `MAR` nor `MNAR`,
  the slices summed to *less* than 100 %. Consistent with `plt_type` and with `plt_group` after the
  `pc()` change, and required by C15 (the old `info > 0` filter is exactly the "NaN becomes invisible"
  pattern).
- `plt_type`'s legend rebuilt its shades from `length(groups)` — the number of feature × group
  *cells*, not of groups — so legend swatches did not match the bars. It now reuses the shades
  actually drawn.

Also fixed in passing: `plt_density`'s `@return` said "before (solid cyan) / after (dashed magenta)"
while the code draws before *dashed* and after *solid*; and the row label `NMV` (no such category) is
now driven off `MV_COLORS`, so the pie's, `plt_group`'s and `plt_type`'s category names cannot drift
apart.

**Validation.** Three tabular fixtures byte-identical; PDF still 5 pages; 8/8 R5 regression tests;
38 new plot assertions (`test_plots.R`) that capture what would actually be *drawn* — slice/colour/
label alignment, `NaN` never reaching a label, matrix shape with a single group, and the `plt_density`
argument contract — each symptom first reproduced against the pre-fix functions (`test_old_symptoms.R`)
so the suite cannot pass vacuously. `lintr` is not installed here, so static checking is
`parse()` plus XML well-formedness only.

### Step 10 (wrapper plumbing and outputs) — implemented

Covered: C4, C19, C20, R12 (`source_local` half), R13.

**One fixture regenerated, deliberately.** `outDmShT1.tabular` is now written unquoted (C20).
§3.1's byte-identity gate therefore applies to steps 2, 5, 6 **and 9**, and to step 10 for SM/VM
only. The regeneration was verified to be a format change and nothing else: the old and new files are
`identical()` once `"` is stripped from every one of the 3076 lines, `identical()` as parsed tables,
max absolute numeric difference **0**, identifier column `identical()`, same 3075 × 24 dimensions. Two
independent re-runs reproduce the new file byte-for-byte, and SM/VM are still byte-identical to their
fixtures. The `min` path was re-run before and after this commit for the same parameters: quoting
only, max difference 0. The only consumer of the fixture is `ImputeMSmv.xml:178`.

**`source_local` swap, verified rather than assumed.** `ImputeMSmv_wrapper.R:26-42` (the `get_arg()`
helper, the `--tool_dir` lookup, the `stop()` guard, and the `base::source()`) is now
`source_local("ImputeMSmv_function.R")`, matching `w4mconcatenate/concatenation_wrapper.R:7-8`, and
`--tool_dir` is gone from the command block. Confirmed on a **relocated** copy of the tool with
`ImputeMSmv_function.R` deleted: the job exits 1 naming
`cannot open file '/tmp/imv_nofile_tool/ImputeMSmv_function.R'` — the script's own directory, not the
job cwd, which is the `get_base_dir()` resolution this finding rests on, and a better message than the
old `stop("tool_dir not provided")`. `tool_dir` no longer appears anywhere in the job log, so
`print(para)` no longer echoes the install path. An audit of all 20 `para$…` keys against the 21 XML
flags now finds **no** key that fails to match a flag exactly.

**C19 was robustness, not a live bug** — worth saying plainly in the reply to Copilot. Pre-fix,
`para$mnar_method` resolved by partial matching, and the `min` log line printed
"imputed using 50 percent of the lowest value…" identically before and after the fix. The hazard is
prospective, exactly as §1 says: adding any sibling key beginning `mnar_method` turns the expression
into `NULL` and `if` into *"argument is of length zero"*. Measured both halves of that.

**One latent bug fixed in passing, by C4's rewrite.** The old
`colnames(x)[colnames(x) == "df"] <- dName` renamed *every* column called `df`, not just the
placeholder — so a sample, feature or metadata field genuinely named `df` came out with two
identically-named columns. Reproduced (`"name" "name" "other"`), and avoided by naming the identifier
column by position.

**Note for the maintainers on `rm(list = ls())` (R13).** Dropped. It is a no-op in a Galaxy job and it
wipes the workspace of anyone who sources the wrapper to debug it. `w4mconcatenate` still has one at
`:1`; if this is adopted as a convention it should be a repo-wide sweep, not left inconsistent.

**Validation.** Gate green (3 fixtures byte-identical against the regenerated set, 5-page PDF); 8/8
R5 tests; 38/38 plot assertions; the 5 negative parameter cases still fail for the intended reasons
(`badcol`, `limnan_empty`, `limnan_range`, `pmin0` each with the specific message, `pmin05` succeeding);
`min` path re-verified. `planemo` is not installed here, so `planemo test` / `planemo lint` /
`planemo shed_lint` have **not** been run — CI must cover them.

*Rig note, since it affects how much the above can be trusted:* the first run through `source_local`
failed with `object 'fake' not found`, which was a defect in the stand-in (its ported helpers looked
each other up in globalenv, where the shim's storage had already been `rm()`'d), not in the tool.
Fixed by giving each ported function the fake namespace as its enclosing environment, mirroring a real
package NAMESPACE. Worth remembering that a shim can fail in ways the real package cannot.

### Step 11 (metadata and docs) — implemented, with two of my own proposals rejected

Covered: C1, C2, C9 (partly), C18, R22, R23, plus two `@param` mismatches from R24 and the
`intesity` typo.

**Three proposals in this report were wrong, and I did not follow them:**

1. **`:cite:` (C2).** §1 proposes rendering the method reference as ``:cite:`doi:...` ``. No tool in
   this repo uses `:cite:` (grep: 0 files), and it is a *Sphinx* role — Galaxy renders help with plain
   docutils, where an unregistered role prints literally. Used the same named hyperlink the file
   already uses for the Software Heritage link, naming the DOI in the link text.
2. **"Align the README version" (C18).** §2 says README `1.0.0` should match the XML
   `1.0.0+galaxy0`. Measured the convention instead: `w4mconcatenate` (README `1.0.0`, XML
   `1.0.0+galaxy1`) and `mzxmlshaper` (README `1.0.0`, XML `1.0.0+galaxy1`) both record the *base*
   version and omit the Galaxy suffix, so `1.0.0` is not a defect. **Left alone.** The `galaxyID`
   half of C18 was real and is fixed: all four READMEs that carry `@galaxyID` match their XML `id`,
   this tool's was the one exception (`imoutemsmv`).
3. **"Use ASCII in `<creator>` attributes" (C9).** `w4mconcatenate/concat.xml:6` has
   `givenName="Mélanie"`, so accented names in `creator` *are* repo practice. Used the accented forms
   already present in this file's own credits block and bibtex (`Cladière`, `Pétéra`).

**C9 deliberately incomplete — Q6 stays open.** The `**Contact:**` paragraph is gone and replaced by
the support-forum line in the repo's existing form (`**Help and support:** <url>`, as in
`mzxmlshaper`), which is the user-facing route the plan asked to keep. `<creator>` gained Cladière and
Pétéra so it names every contributor the credits block names. **No `email=` was added**: I will not
invent an address, and no W4M-wide personal address for these people is attested anywhere in the repo
(the five `email=` instances elsewhere are all other authors'). The only repo-attested W4M address is
`workflow4metabolomics@proton.me` on `tablemerge`'s organization entry; whether that is the right
contact route is a maintainer decision, not mine. `<creator>` also gained a real W4M URL — the old one
pointed at a PDF bitstream in a university repository (`dipot.ulb.ac.be/.../Current_Protocols.pdf`)
rather than the project site. All three URLs were checked live (HTTP 200) before being written.

**Job-log greeting rewritten (R23).** The first line of every job log addressed the user to *MissMDA*,
an unrelated PCA-based MCAR package. `missMDA` appears nowhere else in this branch, so it is a local
slip. Now names ImputeMSmv and describes the MAR/MNAR approach. `Paramters` -> `Parameters`.

**ASCII pass scoped deliberately.** Converted en dashes and `→ ≤ ≠` in help prose and roxygen to
ASCII; the `→` in particular was a silent-corruption risk in the classification rules, where `pNA ≤
LimNan` is a *specification statement*. Kept non-ASCII where it is content, not punctuation: personal
and institutional names (credits, `<creator>`, bibtex), and the `é` in the help's own
`- + * / é` example warning users away from special characters in group names. Post-pass inventory:
`.shed.yml`, wrapper and README are pure ASCII; the XML's 7 remaining lines and function.R's 9 are all
names plus that example. Also added the missing trailing newline to `.shed.yml` (it had none) and
named the method and its DOI in `long_description` (C1), dropping the "Based on the work presented by"
prose that C2 objects to.

**Fixture provenance: not changed, and my own suggestion rejected.** §2 C18 proposed linking the
test-data provenance to "the same Software Heritage origin the help already cites". That link archives
the *original R code*, not a dataset — it cannot identify the fixtures. Writing it into the README
would have manufactured a source. A TODO placeholder naming my own uncertainty is also not something
to ship in documentation, so **the README provenance text is untouched** and this moves to §4 as
**Q7**: which dataset are `inDmShT1/inSmShT1/inVmShT1` derived of? It matters because the tool test
compares all three outputs byte-for-byte. (A second README nit I raised here — "under 1MB" being
wrong because the directory totals 2.9 MB — was mine, and wrong: see Q8, withdrawn after measuring
the files individually.)

**R14 / `digits` not touched.** Listed under step 11 by §3.6, but it is Q1 and changes output
precision and the fixtures. Deferred to the maintainer, consistent with how the P0s were handled.

**Limits of this step's validation.** `planemo` is not installed here, so `planemo lint`,
`planemo xml_schema lint` and `planemo shed_lint` — the actual validators for everything changed in
this commit — **have not been run**. What I could check: XML well-formedness with a parser, `<creator>`
content parsed back out, `.shed.yml` parsed as YAML and confirmed ASCII, RST table widths and rule-line
counts unchanged, R `parse()` on both scripts, byte-identical fixtures, 8/8 R5 and 38/38 plot
assertions, negative cases still failing for the intended reasons. No spell checker (aspell/hunspell)
is installed, so the typo pass covers only the 16 sites found by reading, not a systematic sweep.
