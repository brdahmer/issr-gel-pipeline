# issr-gel-pipeline

Originally built to assist a Master's student project involving population genetics across individuals of a plant species separated geographically in flood-prone (or not) areas using multiple primers. 

Built to somewhat replace some of the subjective steps in the scoring of dominant markers (namely, choosing which bands are "the same band" across lanes and across electrophoresis gels). *Should* automatically measure gel resolution, calculate binning tolerance, pick a stable %tol from a stable plateau, and write a plain `[OK]` / `[CHECK]` checklist for each gel and each primer.

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
| `bin_gel.R`, `combine_gels.R` | R ≥ 4.0. [`ape`](https://cran.r-project.org/package=ape) for bootstrap support on dendrograms (optional) | band table -> loci and 0/1 matrix for one gel (OR several gels of one primer -> one matrix) + population statistics |

Gels must be manually annotated/scored in GelAnalyzer (orignally used with 26.1.1). 
Requires that the user has specified lanes, bands, and ladder sizes (MW), saved as a `.gap` file.

## Quick start / Example

```bash
# STEP 1 - once per gel
python3 extract_gel.py "data/827-g1-superior.gap"
#   -> bands_827-g1-superior.csv, check_827-g1-superior.png, calibration_827-g1-superior.json
#   check_827-g1-superior.jpg: picture of the gel showing bands (yellow), ladder (green), Rf curves (blue), and lane # (red)

# STEP 2 - once per gel
Rscript bin_gel.R "data/bands_827-g1-superior.csv"
#   -> results_827-g1-superior/SUMMARY.txt  (read this first)
#   Optional: --screen, --tol 6.5, --exclude 29, --merge 471+514

# STEP 3 - once per primer (for example, when one primer was run on several gels or populations were divided across gels)
Rscript combine_gels.R "data/samples_827.csv"
#   -> results_827_combined/SUMMARY.txt, combined_matrix.csv
```

See `--help`.

See also **[docs/DOCUMENTATION.md](docs/DOCUMENTATION.md)**.
