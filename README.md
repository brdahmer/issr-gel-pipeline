# issr-gel-pipeline

A small, reproducible pipeline for turning **ISSR gel images scored in GelAnalyzer** into a
**presence/absence (0/1) locus matrix** with population-genetic summaries.
It was built for dominant markers (ISSR, and it works the same way for RAPD and AFLP).

Most scoring of dominant markers still relies on a subjective step: choosing which bands are "the same
band" across lanes and across gels. This pipeline makes that step explicit and testable. It measures the
gel's resolution in pixels, sweeps the binning tolerance, picks a value from a stable plateau, and writes
a plain `[OK]` / `[CHECK]` checklist for each gel and each primer.

```
 GelAnalyzer .gap ──► extract_gel.py ──► bands_<gel>.csv ──► bin_gel.R ──► results_<gel>/
   (one per gel)        STEP 1            + check image        STEP 2        binary_matrix.csv
                                                                                    │
            samples_<primer>.csv  (which lane is which individual / population)    │
                          │                                                         │
                          └──────────────────► combine_gels.R ◄─────────────────────┘
                                                   STEP 3
                                       results_<primer>_combined/
                                       combined_matrix.csv, PhiPT, He, UPGMA, PCoA
```

## Requirements

| Step | Needs |
|---|---|
| `extract_gel.py` | Python ≥ 3.8, **standard library only** (no numpy) |
| `bin_gel.R`, `combine_gels.R` | R ≥ 4.0, base R only. Optional: [`ape`](https://cran.r-project.org/package=ape) for bootstrap support on dendrograms |

Gels are scored in GelAnalyzer (developed with 26.1.1). Each gel needs lanes, bands and ladder sizes
assigned, and is saved as a `.gap` project file.

## Quick start

```bash
# STEP 1 - one call per gel image
python3 extract_gel.py "data/827-g1-superior.gap"
#   -> bands_827-g1-superior.csv, check_827-g1-superior.png, calibration_827-g1-superior.json
#   Open the check image: each cyan curve should follow a row of bands.

# STEP 2 - one call per gel
Rscript bin_gel.R "data/bands_827-g1-superior.csv"
#   -> results_827-g1-superior/SUMMARY.txt  (read this first)
#   Optional: --screen, --tol 6.5, --exclude 29, --merge 471+514

# STEP 3 - one call per primer (only when one primer was run on several gels)
Rscript combine_gels.R "data/samples_827.csv"
#   -> results_827_combined/SUMMARY.txt, combined_matrix.csv
```

Every script prints `--help`. No variables need to be edited inside the scripts.

## Repository layout

```
extract_gel.py            STEP 1  .gap -> band table (Rf, size in bp, intensity)
bin_gel.R                 STEP 2  band table -> loci and 0/1 matrix for one gel
combine_gels.R            STEP 3  several gels of one primer -> one matrix + population statistics
docs/DOCUMENTATION.md     what each script does, its inputs, outputs and options
examples/samples_template.csv   sample-sheet template for step 3
```

See **[docs/DOCUMENTATION.md](docs/DOCUMENTATION.md)** for the technical description.

## Design principles

- **One input file per command.** Defaults are chosen automatically and reported.
- **Nothing is hidden.** Every automatic choice (tolerance, distortion correction, anchor bands) is
  written to the log together with the evidence behind it.
- **Flag, don't silently fix.** Weak lanes, split bands and gel effects are reported for the user to
  decide on. They are never removed automatically.

## License

<!-- Choose a license (e.g. MIT) and add a LICENSE file before publishing. -->

## Citation

<!-- Add how you would like the pipeline to be cited (thesis, paper, or Zenodo DOI). -->
