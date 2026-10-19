# Technical documentation

Three scripts, run in order. Each step reads the previous step's output files, so the steps can be
re-run independently. For example, step 2 can be re-run with another tolerance without redoing step 1.

---

## Step 1 · `extract_gel.py` — from a GelAnalyzer project to a band table

**What it does.** It reads the binary `.gap` file saved by GelAnalyzer: the embedded gel image, the
lanes, the band boundaries and peaks, and the ladder sizes. It then computes, for every band, a
distortion-corrected relative mobility (Rf), a size in base pairs and an intensity.

**How.**

1. **Distortion correction (Rf curves).** Gels are rarely straight: they tilt and "smile". With
   `--curves auto` (default) the script follows each *row* of bands (the same band in neighbouring lanes)
   from lane to lane and fits a line, or a parabola when the curvature is real. These fitted rows serve as
   reference curves of constant Rf, so a band's Rf is interpolated between the curves above and below it.
   When the gel is tilted by more than 4 px end to end, the rows are tracked a second time with the tilt
   removed, so that bands from different rows are not chained together. `--curves drawn` uses the curves
   drawn by hand in GelAnalyzer, and `--curves none` uses Rf = position / lane length.
2. **Size calibration.** Ladder bands give pairs (Rf, bp). Size is interpolated piecewise-linearly in
   log10(bp) between neighbouring ladder bands, so the fit passes exactly through every ladder band.
   Bands outside the ladder range are marked `extrapolated = 1`. `--mw-model exp` fits a single global
   curve `bp = a·exp(−b·Rf) + c` instead.
3. **Densitometry.** For each band it reports peak height and volume above background, computed from
   the embedded image. Lanes whose total volume is under 1/5 of the median are reported as weak (possible
   PCR failure), but they are never removed automatically.

**Outputs** (next to the `.gap`, or in `--out`):

| File | Content |
|---|---|
| `bands_<gel>.csv` | one row per band: `lane, sample, band, ladder, peak_px, rf, mw_bp, extrapolated, raw_volume, rel_volume_pct, …` |
| `check_<gel>.png` | gel crop with sample bands (yellow), ladder bands (green), Rf curves (cyan) and lane numbers (red) |
| `calibration_<gel>.json` | ladder residuals per size, Rf curves used, gel tilt, model used |

`peak_px` is the band position in the *corrected* frame, rounded to whole pixels. Step 2 uses it to
reason in pixels.

**Options:** `--curves auto|drawn|none`, `--out DIR`, `--ladder-lane N`, `--mw-model interp|exp`,
`--auto-tol-px`, `--auto-min-frac`, `--no-tilt-correction`, `--validate GELANALYZER_EXPORT.csv`
(compares Rf and bp against GelAnalyzer's own *Bands Information* table).

---

## Step 2 · `bin_gel.R` — from bands to loci for one gel

**What it does.** It groups bands of similar size into loci ("binning") and writes the 0/1 matrix
(sample × locus). The key decision is the **tolerance**: how different two sizes may be and still count
as the same locus.

**How.**

1. **Clustering.** Band sizes are clustered on log10(bp) with complete linkage, and the tree is cut so
   that no locus spans more than *tol* %. Loci are named `<gel>_<size>bp`.
2. **Resolution in pixels.** Band positions are integer pixels, so one pixel corresponds to a fixed
   percentage of size on that gel. The script reports this value. A tolerance below one pixel cannot
   group anything, so only tolerances of 1 px or more are considered.
3. **Tolerance sweep and plateau.** Tolerances from 0.5 to 15 % are tried. The **plateau** is the
   longest run of consecutive tolerances where (a) the number of loci stays the same and (b) no sample
   has two bands in the same locus. The middle of the plateau is used, unless `--tol` is given.
4. **Quality checks.** These are reported in `SUMMARY.txt`:
   - *merged bands*: one sample with two bands in one locus, which means the tolerance is too large;
   - *physical ceiling*: the occupied lane length divided by the smallest separation between two bands
     of the same lane. The gel cannot resolve more loci than this;
   - *contiguous-pixel cuts*: a locus boundary that falls inside a continuous run of observed
     positions. Such a split is not supported by the data;
   - *linkage agreement*: if complete and single linkage give the same partition, the grouping is unique;
   - *possible split bands*: neighbouring loci ≤ 2 px apart that share no sample, which usually means
     one band class split by gel smiling. They can be joined with `--merge 471+515`;
   - *locus stability*: the fraction of tolerances within ±25 % that keep exactly the same bands in a locus.
5. **Statistics.** For each locus it gives band frequency, expected heterozygosity He (allele frequency
   q = √(1 − band frequency), as in GenAlEx), Shannon's I, and PIC = 2f(1 − f). It also computes a
   Jaccard distance matrix and a UPGMA dendrogram with its cophenetic correlation. If `ape` is installed,
   the dendrogram gets bootstrap support (1000 resamples of loci).

**Outputs** (`results_<gel>/`):

| File | Content |
|---|---|
| `SUMMARY.txt` | plain-language checklist — read first |
| `binary_matrix.csv` / `.txt` | the 0/1 matrix (main result) |
| `loci_bins.csv` | per locus: mean/min/max size, span, gap to neighbours, conflicts |
| `locus_statistics.csv` | per locus: frequency, He, I, PIC, stability |
| `primer_summary.csv` | one-line summary (used by step 3) |
| `tolerance_sweep.csv` | the full sweep table |
| `bands_annotated.csv` | input bands with their assigned locus |
| `jaccard_distance.csv`, `dendrogram_*.pdf`, `upgma_tree_*.nwk` | sample clustering |
| `diagnostics_*.pdf`, `matrix_*.pdf` | sweep plots, band map, matrix heat map |
| `analysis_log.txt` | parameters, session info |

**Options:** `--screen` (print the sweep table only), `--tol X`, `--exclude 9,10`, `--merge A+B[,C+D]`,
`--keep-outside-ladder`, `--out DIR`.

---

## Step 3 · `combine_gels.R` — several gels of one primer into one matrix

**What it does.** When the samples of one primer were run on several gels (or on the upper and lower
tiers of one gel), it matches the loci of the different gels by size. It then builds one combined
matrix, tests whether the gel itself left a signature in the data, and computes population statistics.

**Input.** A sample sheet (`examples/samples_template.csv`) with columns `tier, lane, individual,
population, notes`. `tier` is the gel name used in steps 1–2. `lane` is the **red** lane number on
`check_<tier>.png`. A row is used only when `individual` is filled in and `notes` is empty, so write
`ladder`, `failed PCR`, etc. in `notes` to leave a lane out. The script looks for `results_<tier>/` and
`calibration_<tier>.json` next to the sheet, or in `--dir`.

**How.**

1. **Locus matching.** Loci from different tiers are clustered on log(bp) with complete linkage, under
   one constraint: two loci of the same tier are never merged, because step 2 already decided they differ.
2. **Anchor bands and tolerance.** Loci present in ≥ 50 % of the samples of two or more tiers and within
   5 % of each other are treated as the same band. Their largest size difference estimates the
   calibration error between gels, and the matching tolerance is set to 1.5 × that error (minimum 1 %).
   Without anchors, the middle of the longest stable plateau is used instead. When no anchors are found,
   the script also checks whether one constant size shift would align the frequent bands of two gels.
   That correction can be applied deliberately with `--shift tier=−9.4`.
3. **Replicates.** An individual listed on two tiers enters the matrix once. The percentage of loci that
   agree between its two profiles is reported.
4. **Gel-effect test.** When a population is spread over two or more tiers, PhiPT (AMOVA) is computed
   between tiers within that population, along with Fisher's exact test per locus. A significant result
   means either a gel artefact or real substructure. The summary also warns when tiers and populations
   are fully confounded.
5. **Population statistics.** Per population it reports % polymorphic loci, He, Shannon's I and
   bands per individual. It also computes PhiPT between populations (999 permutations), all repeated at
   the neighbouring tolerances as a sensitivity check. It draws a UPGMA/Jaccard dendrogram (with
   bootstrap support when `ape` is installed) and a PCoA, coloured by population and shaped by tier.

**Outputs** (`results_<primer>_combined/`): `SUMMARY.txt`, `combined_matrix.csv` / `.txt` (individual,
population, loci), `combined_loci.csv` (which tier loci were merged, frequency per population and tier),
`population_diversity.csv`, `gel_effect_test.csv`, `tolerance_sensitivity.csv`, `dendrogram_*.pdf/png`,
`upgma_tree_*.nwk`, `pcoa_*.pdf/png`, `combination_log.txt`.

**Options:** `--tol X`, `--dir DIR`, `--out DIR`, `--name PRIMER`, `--permutations N`,
`--shift tier=percent[,tier=percent]`.

---

## Notes and limitations

- Band *calls* come from GelAnalyzer and the user. The pipeline measures and bins them, but it does
  not detect bands itself.
- Without replicate samples across gels, the genotyping error rate cannot be estimated. Cross-gel
  locus matching then depends only on each gel's ladder calibration.
- Statistics assume dominant markers and Hardy–Weinberg equilibrium for allele frequencies (the usual
  convention for ISSR, as in GenAlEx).
