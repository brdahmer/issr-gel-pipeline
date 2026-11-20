# issr-gel-pipeline

Originally built for a Master's project on the population genetics of a plant species, comparing individuals from flood-prone and non-flooded areas with several ISSR primers.

Scoring dominant markers involves some subjective calls, mainly deciding which bands are "the same band" across lanes and across gels. The scripts try to make those calls from the gel itself: they measure the gel's resolution, try a range of binning tolerances, pick one from the middle of a stable range, and write an `[OK]` / `[CHECK]` checklist for each gel and each primer.

```
 GelAnalyzer .gap file ──► extract_gel.py ──► bands_<gel>.csv ──► bin_gel.R ──► results_<gel>/
      (one per gel)           STEP 1            + check image       STEP 2       binary_matrix.csv
                                                                                    │
            samples_<primer>.csv  (lane to individual/population correspondence).   │
                          │                                                         │
                          └──────────────────► combine_gels.R ◄─────────────────────┘
                                                   STEP 3
                                       results_<primer>_combined/
                                       combined_matrix.csv, PhiPT, He, UPGMA, PCoA
```

## Requirements

| Step | Needs | Input -> Output |
|---|---|---|
| `extract_gel.py` | Python ≥ 3.8 | .gap -> band table (Rf, size in bp, intensity) |
| `bin_gel.R`, `combine_gels.R` | R ≥ 4.0; optionally [`ape`](https://cran.r-project.org/package=ape) for bootstrap support on the dendrograms | band table -> loci and 0/1 matrix for one gel (or several gels of one primer -> one matrix), plus population statistics |

Gels must first be scored by hand in GelAnalyzer (developed with version 26.1.1): mark the lanes, bands and ladder sizes, and save the `.gap` file.

## Quick start / Example

```bash
# STEP 1 - once per gel
python3 extract_gel.py "examples/827-g1-superior.gap"
#   -> bands_827-g1-superior.csv, check_827-g1-superior.png, calibration_827-g1-superior.json
#   check_827-g1-superior.png: picture of the gel showing bands (yellow), ladder (green), Rf curves (blue), and lane # (red)

# STEP 2 - once per gel
Rscript bin_gel.R "examples/bands_827-g1-superior.csv"
#   -> results_827-g1-superior/SUMMARY.txt  (read this first)
#   Optional: --screen, --tol 6.5, --exclude 29, --merge 471+514

# STEP 3 - once per primer (for example, when one primer was run on several gels or populations were divided across gels)
Rscript combine_gels.R "data/samples_827.csv"
#   -> results_827_combined/SUMMARY.txt, combined_matrix.csv
```

Each script lists its options with `--help`. [DOCUMENTATION.md](DOCUMENTATION.md) explains how each step works and what it writes.
