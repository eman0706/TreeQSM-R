# TreeQSM-R

An R port of **TREEQSM** (Raumonen et al., 2013; original MATLAB code at https://github.com/InverseTampere/TreeQSM), the toolbox that builds quantitative structure models (QSMs) of trees from point clouds,
together with the code used to validate it on under-canopy UAV structure-from-motion point clouds of a Mediterranean Aleppo pine stand.

> **Status: research code.** This is not an official release of TreeQSM and is not endorsed by its authors. The state tagged `v0.2.0-submission` is the one used for the results of the accompanying paper (see *Citation*).

## What has been validated, and what has not

* **Stem diameters (DBHcyl at breast height): compared with the original TreeQSM 2.4.1** run in GNU Octave 10.1.0 on the same 41 extracted trees (two point clouds) and the same eight parameter combinations. Per-tree differences between the port and the original: median 0.02 cm (cloud 15, 21 trees) and 0.01 cm (cloud 11, 20 trees), largest 0.90 and 0.27 cm, all trees within 1 cm; the differences are of the order of the variability between repeated runs of the original itself (random numbers cannot be seeded identically in the two programs).
* **Field validation:** median absolute error of DBHcyl against caliper diameters 3.4% and 3.1% (original in Octave: 3.4% and 2.9%).
* **Not validated:** tree height, number of cylinders, branch structure, volumes of branches and total tree volume. The number of segments produced by the two segmentation passes differs from the original by about a factor of two in the diagnosed trees. **Do not use the branch attributes of this port without your own comparison with the original.**
* **Known deviations** from TreeQSM 2.4.1 are listed in Table 10 of the paper (random numbers, fixed orthonormal basis, `define_trunk`, branch identifiers and order, radius cap, single-precision handling, gap filling not ported).
* Tested with Octave, not with MATLAB.

## Contents

| folder | content |
|---|---|
| `R/` | the library (files that define the functions reached by `treeqsm()`), the faithful versions `define_base_forb_fedele` (in `prova_tree_sets_fedele.R`) and `least_squares_cylinder_fedele.R`, the clean-session loader `sessione_pulita.R`, and the analysis scripts of the paper |
| `octave/` | driver that runs the **unmodified** original TreeQSM in GNU Octave on the extracted trees, plus diagnostic scripts used to compare stages of the two programs |
| `python/` | auxiliary analyses in Python (solver comparison, direct estimator, call graph); paths must be adapted |
| `docs/` | `GLOSSARY.md` (English guide to the scripts), annotated call tree, inventory of the functions, manifest of the validated library (function, defining file, md5), `sessionInfo.txt` |

The point clouds, the field measurements and the per-tree results of other authors are **not** included: they belong to R. Greco and P. Tartarino and are available from them on reasonable request.

## Language note

The names of the original TREEQSM functions are in English, as in MATLAB. The analysis and validation scripts, their function names and the messages printed to the console are in **Italian**, and were left unchanged on purpose, because the state tagged `v0.2.0-submission` is the one validated in the paper. `docs/GLOSSARY.md` lists every script, its main function and what it does, in English.

## Requirements

R 4.6.1 on Windows 10 (Intel Core i9-7900X, 64 GB RAM). The QSM functions and the diameter estimators use only base R. The packages `lidR` (4.3.2) and `dbscan` (1.2.6) are needed only to read LAS files and to isolate the trees. The code uses `%||%` (base R >= 4.4.0). GNU Octave 10.1.0 with the `statistics` package (Octave Forge) is needed only for the comparison with the original.

## Load the final configuration from a clean R session

```r
# new R session, nothing restored (check that ls() is empty)
source("R/sessione_pulita.R")
carica_libreria_pulita(base = "R")   # fedele = TRUE: faithful base of the second pass and Gauss-Newton solver of the original
```

This loads an explicit list of files, activates `define_base_forb_fedele` and `least_squares_cylinder_fedele`, and writes `manifesto_sessione_pulita.csv`; `confronta_manifesto()` compares it with `docs/manifesto_validato.csv`.

## Reproduce the numbers of the paper (needs the data)

```r
# 1. models (8 per tree, 3 seeds), per cloud
raccolta_modelli_quote(<folder of the extracted trees>, offsets = 0:2, file_out = "modelli_quote_fedele.rds")
# 2. summaries, paired tests, comparison with the original
analisi_fedele(c15, c11)
tabella7_fedele(R)
confronta_octave_fedele(c15, c(<octave csv rip2>, <octave csv rip1>), "cloud 15")
# 3. clean-session reproducibility
riproduci_fedele(base = "R")
```

Running times on the machine above (sequential, one R process): 39 min (cloud 15, 21 trees) and 62 min (cloud 11, 20 trees) for 8 models x 3 seeds per tree.

## Comparison with the original in Octave

See `octave/GUIDA_OCTAVE.md`. The driver `esegui_treeqsm_octave.m` runs the original source files unchanged; the only addition is a stub for `alphaShape` (not available in Octave), which affects only the crown measures based on alpha shapes (not used here, returned as NaN).

## Use of an AI assistant

The code, the tests and a first draft of the paper were developed with an AI assistant (Claude, Anthropic): translation of the MATLAB functions into R, corrected and faithful versions of some functions, diagnostic tools and test scripts. All R and Octave code was run by the authors; the assistant never executed it. See Section 3.6 and the declaration of the paper.

## License and attribution

Released under the **GNU General Public License, version 3 or later**, as the MATLAB sources from which it derives (each original file states "This file is part of TREEQSM" and the GPL v3-or-later terms). See `LICENSE` and `NOTICE.md`.

## Citation

See `CITATION.cff`. If you use this code, please also cite Raumonen et al. (2013), *Remote Sensing* 5(2), 491-520, doi:10.3390/rs5020491, and Greco et al. (2023), *Int. J. Appl. Earth Obs. Geoinf.* 122, 103426, doi:10.1016/j.jag.2023.103426.
