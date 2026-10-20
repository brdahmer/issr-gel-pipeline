## Step 1 · `extract_gel.py`

Reads the `.gap` file from GelAnalyzer after manual scoring/annotation. Typically includes/needs an embedded gel image, lane and band selections, well (0) and bottom (1) Rf curves, ladder info and attributed MW, then calculates relative mobility (Rf curves), size in bp and intensity (densitometry).

**How.**

1. **Distortion correction (Rf curves):** Corrects natural tilt and "smile" events in electrophoresis gels. With
   `--curves auto` (default), the script follows each *row* of bands (or similar height band in neighbouring lanes)
   from lane to lane and fits either a line (1st degree or 2nd degree when the curvature is real). Fitted rows serve as
   reference curves of constant Rf.
   When the gel is tilted by more than 4 px end to end, the rows are tracked a second time with the tilt
   removed, so that bands from different rows are not chained together.
     - `--curves drawn` uses the curves drawn by hand in GelAnalyzer
     - `--curves none` uses Rf = position / lane length.
   
2. **Size calibration:** Estimates bp using the physical ladder lane by pairing the known bp dimension (from ladder) with their calculated Rf.
   By default (`--mw-model interp`), the script connects these points using local piecewise-linear interpolation using log10(bp)
   between neighboring ladder bands. Any bands migrating beyond the ladder range are marked `extrapolated = 1`.
     - `--mw-model exp` fits a single global curve `bp = a·exp(−b·Rf) + c` instead if the gel fits an exponential curve.
   
3. **Densitometry:** Analyzes the embedded gel picture to calculate peak height and background-subtracted volume for every lane. If a lane volume falls below 20% of the median lane intensity, the lane is reported as weak (*possible PCR failure*) but are not automatically removed.

**Outputs** (files are saved directly alongisde the `.gap` file or in a folder specified with `--out`):

| File | Content |
|---|---|
| `bands_<gel>.csv` | primary output in which every row reports a single band's metric: `lane, sample, band, ladder, peak_px, rf, mw_bp, extrapolated, raw_volume, rel_volume_pct, etc.` |
| `check_<gel>.png` | gel crop with sample bands (yellow), ladder bands (green), Rf curves (blue) and lane numbers (red) |
| `calibration_<gel>.json` | ladder residuals per size, Rf curves used, gel tilt angles, model (regression constants) used |

`peak_px` is the band position in the *corrected* frame, rounded to whole pixels. Step 2 uses it to
reason in pixels.

**Command-Line Arguments**
```
Core execution syntax:
python extract_gel.py "gel-data.gap" --curves auto --mw-model interp

Custom flags:
--curves [auto|drawn|none]    # Pick geometric correction behavior
--mw-model [interp|exp]       # Pick size interpolation vs. exponential curve fitting
--ladder-lane N               # Specify a custom lane as ladder
--out DIR                     # Set custom output destination
--no-tilt-correction          # Force-disable automated tilt adjustment
--validate EXPORT.csv         # Compares Rf and bp against GelAnalyzer's own Bands Information table
--auto-tol-px / --auto-min-frac # Finetune row tracking thresholds
```
---

## Step 2 · `bin_gel.R`

Classified clean and sized bands into genetic markers (loci). Groups bands of similar size across all samples into loci ("binning") and writes the binary presence/absence (1/0) matrix for downstream population genetics analysis.

**How.**

1. **Clustering:** Runs one-dimensional complete linkage hierarchical clustering on the log10(bp) size of bands and cuts the tree so that no locus spans wider than a set % tolerance (`tol`). Tolerance values <1 px are ignored.

2. **Tolerance scan and plateau:** Scans across `tol` values from 0.5% to 15%, and defines the stable plateau as the longest run of consecutive `tol` values where (a) the number of loci remains stable and (b) no single individual sample has two of its own bands compressed into the same locus bin. The middle of the plateau is automatically employed unless `--tol` is used.

3. **Quality checks:** Writes a report in `SUMMARY.txt`;
   - *merged bands*: Warns if `tol` is too high, causing two distinct bands from the same individual to fuse into one locus.
   - *physical ceiling*: Computes a physical limit based on the gel resolution; if the locus count exceeds the limit, warns the user about over-interpreting noise.
   - *contiguous-pixel cuts:* Flags instances where a locus boundary was forced right through a continuous cluster of pixels, which is not supported by data. 
   - *linkage convergence*: Cross-checks complete linkage vs. single linkage (maximum vs. minimum distance). If both are identical, the locus boundaries are hopefully mathematically robust.
   - *split band detection*: Flags neighboring loci separated by <2 px that share no overlapping sample, which usually means an uncorrected gel smiling, and can be joined with `--merge A+B[,C+D]`;
   - *locus stability profile*: Tests every locus by fluctuating the optimal `tol` by +-25%.
   
4. **Baseline population metrics.**
By threating markers as a dominant genetic system, it calculates band frequencies, Expected Heterozygosity (He), Shannon's Information Index (I), and Polymorphism Information Content (PIC). Then computes a global Jaccard distance matrix and builds an UPGMA dendrogram. If the R package `ape` is available, it automatically appends confidence values (1,000 locus bootstrap resamples).

**Outputs** (outputs are written inside a dedicated folder `results_<gel>/`):

| File | Content |
|---|---|
| `SUMMARY.txt` | checklist summarizing biological and physical quality of the run |
| `binary_matrix.csv` / `.txt` | binary presence/absence matrix |
| `loci_bins.csv` | locus structural data, mean/min/max size, span (%width), separation gap to nearest neighbours, sample conflicts |
| `locus_statistics.csv` | locus genetic data, band frequency, He, I, PIC, stability |
| `primer_summary.csv` | single-row summary table containing global quality for Step 3 if necessary |
| `tolerance_sweep.csv` | full tolerance sweep |
| `bands_annotated.csv` | original input dataset appended with the final locus bin assignments |
| `jaccard_distance.csv`, `dendrogram_*.pdf`, `upgma_tree_*.nwk` | sample clustering, distance matrices, UPGMA tree plots, and Newick tree files |
| `diagnostics_*.pdf`, `matrix_*.pdf` | multi-panel visual plots containing the tolerance sweep, spatial band map, binary matrix heat map |
| `analysis_log.txt` | execution log, parameters, session info |

**Command-Line Arguments**
```
Core execution syntax
Rscript bin_gel.R "bands_gel_i5.csv"

# Custom tweaks
--screen               # Run a quick trial; prints the sweep table only then exits
--tol 6.5              # Bypass auto-optimization and enforce a fixed percentage tolerance
--exclude 9,10         # Completely drop specific low-quality lanes from the matrix
--merge 471+515        # Force-fuse artifactual split bins into a single locus
--keep-outside-ladder  # Prevent the script from discarding extrapolated bands
--out DIR              # Redirect all output files to a custom directory 
```
---

## Step 3 · `combine_gels.R`

When samples for (for example) a single primer span multiple different gels, or are split between the upper and lower tiers of the same gel, sizes can drift. The script attempts to match up separated loci by size to build a single, master binary matrix, then executes the binning and statistical checks.

**Input.** 

Requires a sample sheet (`examples/samples_template.csv`) with the following column headers: 
- `tier`: exact gel name used in steps 1 and 2
- `lane`: must correspond to the red lane number on `check_<gel>.png`. 
- `individual`: unique ID of each sample. 
- `population`: population or collection group name     
- `notes`: used to filter out lanes/noise. If anything is typed in the `notes` column (like _ladder_, _failed_PCR_, or _blurry_, for example), the script skips the lane. Only the rows in which `individual` is filled in and `notes` is blank are processed.

The script expects to find `results_<gel>/` folders and `calibration_<gel>.json` files sitting in the working directory alongside the sample spreadsheet (or use `--dir`).

**How.**

1. **Locus matching:** Runs the complete linkage hierarchical clustering algorithm on the log10(bp) sizes and clusters loci from different gels. As Step 2 should have already determined what bands are distinct, the script does not merge two loci of the same gel.
   
2. **Calibration Error (Anchor Bands) and tolerance:** To determine how much band sizes drift between different gel runs, the script searches for _anchor bands_, _i.e._ loci present in >50% of the samples in two or more gels and within 5% of each other.
     - The script finds the single largest size difference among anchors to estimate cross-gel calibration error, then automatically sets cross-gel matching tolerance to 1.5 × that maximum error (baseline minimum of 1%).
     - If gels do NOT share common high-frequency bands, the script scans for a stable plateau across a range of tolerances (as before) and picks the middle. It will also calculate whether a fixed mathematical shift would line up frequent bands. If a systemic drift is noticeable, the user can force an adjustment using `--shift tier=-9.4`, for example.

3. **Replicate handling:** If one individual is ran on both gels as a quality check, it is compressed into a single entry in the final matrix, then calculates the %agreement between the two gel runs.

4. **Gel-effect test:** If a population is spread over two or more gels, the script extracts those specific groups and runs an AMOVA to calculate PhiPT between the gels alongside a locus-by-locus with Fisher's test.
     - If this returns a highly significant p-value, warns the user that it likely reflects physical gel artifact.
     - Also warns the user if populations and gels are completely confounded, _e.g.,_ Population 1 was ran on Gel A and Population 2 was ran on Gel B, etc., and separating biological from mechanical differences is relatively impossible.

5. **Population metrics:** As before from the combined matrix, then computes overall PhiPT across populations using 999 random permutations and tests nearby `tol` values to ensure conclusions don't fall apart if the clustering threshold shifts. Finally, plots an UPGMA/Jaccard Dendrogram as before and produces a PCoA scatter plot.

**Outputs** (written to `results_<primer>_combined/`):
| File | Content |
|---|---|
| `SUMMARY.txt` | |
| `combined_matrix.csv` / `.txt` | unified presence/absence matrix organized by individual and population |
| `combined_loci.csv` | mapping showing which individual gel bins were merged, including their specific locus frequencies per population and gel |
| `population_diversity.csv` | summary table of standard genetic diversity indicators |
| `gel_effect_test.csv` | do the individual physical gel runs skew data? |
| `tolerance_sensitivity.csv` | total locus counts and population PhiPT scores across alternative thresholds |
| `dendrogram_*.pdf/png` & `upgma_tree_*.nwk` | visual tree diagrams and raw Newick tree files mapping individual sample relationships |
| `pcoa_*.pdf/png` | scatter plot designed to visually isolate population grouping from gel-run cluster patterns |
| `combination_log.txt` | log of the run |


**Command-Line Arguments**

**Command-Line Arguments**
```
Basic cross-gel merge execution
Rscript combine_gels.R "my_sample_sheet.csv"

# Advanced parameters
--tol X                      # Force a specific cross-gel size tolerance and bypass auto-anchoring
--dir DIR                    # Location of Step 2 folders and JSON files
--out DIR                    # Redirect all outputs to a specific directory
--name PRIMER                # Define a specific name tag for your outputs
--permutations N             # Adjust the total number of random shuffles for AMOVA significance
--shift tier=percent[,tier=percent] # Intentionally shift specific gel sizes by a percentage to correct drift
```

---

## Notes and limitations

- NOT an automated peak detector, as it relies entirely on the initial visual bands and manual annotations done in GelAnalyzer. If a true biological band was missed or a visual artifact was misidentified during that initial phase, the pipeline will bin what it is given.
- Cross-gel locus alignment is highly dependent on the precision of each gel's unique ladder calibration curve. _Without running explicit replicate samples across different gels, the pipeline CANNOT reliably calculate a statistical genotyping error rate_. Please add a few duplicate individuals across separate runs as (by standards) highly recommended to enable absolute error parameters calculation.
- All downstream diversity metrics, heterozygosity indices, and population statistics are computed under the usual convention for ISSR for dominant molecular markers like in GenAlEx.
