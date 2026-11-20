## Step 1 · `extract_gel.py`

Reads a `.gap` file saved by GelAnalyzer after you have marked the lanes, bands and ladder sizes. The file holds the gel image, the lanes and bands, the well (Rf 0) and bottom (Rf 1) curves and the ladder sizes. The script works out each band's relative mobility (Rf), its size in bp and its intensity.

### How it works

1. **Tilt and smile correction.** With `--curves auto` (the default), the script follows each row of bands (the same band at a similar height in neighbouring lanes) across the gel. It fits a straight line to each row, or a parabola when the row is clearly curved, and uses the fitted rows as curves of constant Rf. If the gel is tilted by more than 4 px from one side to the other, it tracks the rows a second time with the tilt removed, so bands from different rows are not chained together.
   - `--curves drawn` uses the curves drawn by hand in GelAnalyzer.
   - `--curves none` uses Rf = position / lane length.

2. **Size calibration.** Each ladder size is paired with its Rf. By default (`--mw-model interp`) the size of a band is interpolated in log10(bp) between the two ladder bands on either side of it. Bands beyond the ladder range get `extrapolated = 1`. `--mw-model exp` fits one curve, `bp = a·exp(−b·Rf) + c`, to the whole ladder instead.
   - Ladder sizes must decrease down the gel. If they don't (usually a size typed on the wrong band), the script stops and names the sizes.
   - With `--ladder-lane N`, only lane N is used to calibrate.

3. **Intensity.** The script reads the gel image and measures each band's peak height and its volume above the lane background. A lane whose total volume is below 20 % of the median lane is reported as weak (possibly a failed PCR); a lane with no band marked counts as 0. Weak lanes are reported, not removed.

If lanes have different lengths, the script warns you, because step 2 compares pixel positions across lanes.

### Outputs

Files are written next to the `.gap` file, or to the folder given with `--out`.

| File | Content |
|---|---|
| `bands_<gel>.csv` | the main output, one row per band: `lane, sample, band, ladder, peak_px, rf, mw_bp, extrapolated, raw_volume, rel_volume_pct`, etc. |
| `check_<gel>.png` | the gel with sample bands (yellow), ladder bands (green), Rf curves (blue) and lane numbers (red) |
| `calibration_<gel>.json` | ladder residuals per size, the Rf curves used, the gel tilt and the size model |

`peak_px` is the band position after the tilt and smile correction, rounded to whole pixels. Step 2 uses it for its pixel-based checks.

### Command-line options
```bash
python3 extract_gel.py "gel-data.gap"
```

```
--curves auto|drawn|none        tilt and smile correction (default: auto)
--mw-model interp|exp           size calibration: interpolation (default) or exponential fit
--ladder-lane N                 use lane N as the ladder (default: every lane with sizes assigned)
--out DIR                       output folder (default: the folder of the .gap)
--no-tilt-correction            do not remove the gel tilt when tracking rows of bands
--validate EXPORT.csv           compare Rf and bp with GelAnalyzer's "Bands Information" export
--auto-tol-px, --auto-min-frac  row-tracking settings (defaults 3 px and 0.5 of the lanes)
```

---

## Step 2 · `bin_gel.R`

Groups bands of similar size across samples into loci ("binning") and writes the 0/1 presence/absence matrix for one gel.

### How it works

1. **Clustering.** Bands are clustered by log10(bp) with complete linkage, and the tree is cut so that no locus is wider than a tolerance (`tol`, in % of size).

2. **Choosing the tolerance.** The script tries every `tol` from 0.5 % to 15 % in steps of 0.25 %, and looks for the longest run of consecutive values where:
   - (a) the number of loci stays the same,
   - (b) no sample has two of its own bands in one locus,
   - (c) `tol` is no wider than the detection floor, which is the closest pair of bands within one lane, converted from px to %,
   - (d) the number of loci does not exceed what the gel can resolve.

   A run needs at least 3 values (0.5 %) to count. The middle of the run is used unless you set `--tol`. If no run meets (c) and (d), the longest run meeting (a) and (b) is used, with a warning. The script prints the lane and bands that set the detection floor: a band marked twice in GelAnalyzer would lower the floor for the whole gel.

   One pixel is worth a different % of size at different heights on the gel (more near the top and bottom), so the script measures it between each pair of neighbouring ladder bands. Rule (c) uses the finest value found where there are bands. Tolerances narrower than 1 px at the coarsest value are skipped, since they cannot group anything.

3. **Checks.** `SUMMARY.txt` lists each check as `[OK]` or `[CHECK]`:
   - *Merged bands:* `tol` is so wide that two bands from the same sample fall in one locus.
   - *Locus ceiling:* more loci than the gel can resolve.
   - *Contiguous-pixel cuts:* a locus boundary falls inside a continuous run of band positions, where the data show no gap.
   - *Linkage agreement:* complete and single linkage give the same loci. If they do, the grouping does not depend on the method.
   - *Split bands:* two neighbouring loci less than 2 px apart that share no sample. This usually means gel smiling that was not fully corrected; join them with `--merge A+B[,C+D]`.
   - *Locus stability:* whether each locus keeps the same bands when `tol` changes by ±25 %.
   - *Weak lanes:* lanes with very little signal, where missing bands may be failed PCRs.

4. **Diversity.** Treating the bands as dominant markers, the script reports band frequency, expected heterozygosity (He), Shannon's index (I) and polymorphism information content (PIC). It builds a Jaccard distance matrix and a UPGMA dendrogram. If the R package `ape` is installed, it adds bootstrap support from 1,000 resamples of the loci.

### Outputs

Written to `results_<gel>/`.

| File | Content |
|---|---|
| `SUMMARY.txt` | the checklist; read this first |
| `binary_matrix.csv` / `.txt` | 0/1 presence/absence matrix |
| `loci_bins.csv` | per locus: mean, min and max size, width (%), gap to the next smaller locus (% and px), samples with two bands in it |
| `locus_statistics.csv` | per locus: band frequency, He, I, PIC, stability |
| `primer_summary.csv` | one-row summary of the gel, used by step 3 |
| `tolerance_sweep.csv` | results at every tolerance tried |
| `bands_annotated.csv` | the band table with the locus of each band |
| `jaccard_distance.csv`, `dendrogram_*.pdf`, `upgma_tree_*.nwk` | distances between samples, dendrograms and the tree in Newick format |
| `diagnostics_*.pdf`, `matrix_*.pdf` | plots of the tolerance sweep, the bands by size and the 0/1 matrix |
| `analysis_log.txt` | parameters, results and R session details |

### Command-line options
```bash
Rscript bin_gel.R "bands_gel_i5.csv"
```

```
--screen               print the loci at each tolerance and the suggested one, then stop
--tol 6.5              use this tolerance (%) instead of the automatic choice
--exclude 9,10         leave these lanes out (e.g. failed PCR)
--merge 471+515        join these loci (sizes in bp) into one; several groups: 471+515,1333+1399
--keep-outside-ladder  keep bands outside the ladder range (dropped by default)
--out DIR              output folder (default: results_<gel>/ next to the band file)
```

---

## Step 3 · `combine_gels.R`

Combines the gels of one primer into one matrix, for example when the samples were run on several gels, or on the upper and lower rows of one gel. Band sizes drift a little between gels, so the script matches loci across gels by size, then runs the checks and statistics on the combined matrix.

### Input

A sample sheet (see `examples/samples_template.csv`) with these columns:
- `tier`: the gel name used in steps 1 and 2.
- `lane`: the red lane number on `check_<gel>.png`.
- `individual`: a unique ID for each sample.
- `population`: the population or collection group.
- `notes`: anything typed here (e.g. _ladder_, _failed PCR_, _blurry_) leaves the lane out. Only rows with `individual` filled in and `notes` empty are used.

The script looks for the `results_<gel>/` folders and `calibration_<gel>.json` files next to the sample sheet, or in the folder given with `--dir`. It stops if one lane is assigned to two individuals.

### How it works

1. **Matching loci.** Loci from different gels are clustered by log(bp) with complete linkage. Two loci from the same gel are never merged, since step 2 already decided they are different.

2. **Matching tolerance.** The script looks for anchor bands: loci present in at least 50 % of the samples on two or more gels, and within 5 % of each other in size. The largest size difference among the anchors is taken as the calibration error between gels, and the matching tolerance is set to 1.5 × that error (at least 1 %).
   - Without anchors, it tries tolerances from 0.5 % to 10 % and uses the middle of the longest stable run.
   - It also checks whether shifting one gel's sizes by a fixed % would line up its frequent bands with another gel's. If so, you can apply the shift with `--shift tier=-9.4`.

3. **Replicates.** An individual run on more than one gel appears once in the final matrix. The log reports the % of loci scored the same in each copy.

4. **Gel effect.** If a population was run on two or more gels, the script compares those gels within the population: PhiPT from an AMOVA, and Fisher's exact test per locus (p-values adjusted across loci with Benjamini-Hochberg; loci are flagged at q < 0.05). A clear difference between gels is likely a gel artefact. If every gel holds a single population, gel and population cannot be told apart, and the summary says so.

5. **Population statistics.** On the combined matrix: diversity per population, PhiPT between populations (999 permutations by default) and the same statistics at neighbouring tolerances, to check that the conclusions hold. It also draws a UPGMA/Jaccard dendrogram and a PCoA plot.

### Outputs

Written to `results_<primer>_combined/`.

| File | Content |
|---|---|
| `SUMMARY.txt` | the checklist; read this first |
| `combined_matrix.csv` / `.txt` | 0/1 matrix by individual and population |
| `combined_loci.csv` | which gel loci were merged into each combined locus, with band frequencies per population and per gel |
| `population_diversity.csv` | diversity statistics per population |
| `gel_effect_test.csv` | PhiPT and differing loci between gels, within each population |
| `tolerance_sensitivity.csv` | number of loci and PhiPT at neighbouring tolerances |
| `dendrogram_*.pdf/png`, `upgma_tree_*.nwk` | dendrograms and the tree in Newick format |
| `pcoa_*.pdf/png` | PCoA plot, coloured by population, with one symbol per gel |
| `combination_log.txt` | log of the run |

### Command-line options
```bash
Rscript combine_gels.R "samples_827.csv"
```

```
--tol 4                  matching tolerance (%) between gels, instead of the automatic choice
--dir DIR                folder with the step 1 and 2 outputs (default: the folder of the sheet)
--out DIR                output folder (default: results_<primer>_combined/)
--name PRIMER            name used in output files (default: from the sheet's file name)
--permutations N         permutations for the PhiPT tests (default 999)
--shift tier=-9.4[,...]  shift one gel's sizes by a % before matching, to correct calibration drift
```

---

## Notes and limitations

- The pipeline does not detect bands. It bins the bands you marked in GelAnalyzer, so a missed band or a marked artefact goes through unchanged.
- Matching loci across gels depends on how good each gel's ladder calibration is. Without the same individuals run on more than one gel, the genotyping error rate cannot be estimated. Running a few individuals on two gels is strongly recommended.
- Diversity statistics follow the usual convention for dominant markers such as ISSR, as in GenAlEx. He assumes Hardy-Weinberg equilibrium; for selfing or clonal species, rely more on % polymorphic loci, band frequencies, Shannon's I and PhiPT.
