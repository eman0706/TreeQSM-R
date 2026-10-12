# Glossary: file and function names in Italian

The code was developed by Italian-speaking authors. The **names of the original TREEQSM functions** (`treeqsm`, `cover_sets`, `cylinders`, ...) are kept as in MATLAB and are in English. The **scripts of the analysis and the validation**, their function names and the messages printed to the console are in Italian. They were left unchanged on purpose: the state tagged `v0.2.0-submission` is the one validated and cited in the paper (names in its Appendix and Table 10), and the manifest `docs/manifesto_validato.csv` refers to these names. An English-named version may follow in a later release.

Frequent words in the identifiers: *albero* = tree, *alberi* = trees, *nuvola* = point cloud, *quota* = height above ground, *fedele* = faithful (to the original), *pulita* = clean, *prova/prove* = test(s), *diagnosi* = diagnosis, *porting* = port, *campo* = field, *fusto* = stem, *radice* = root (here: the root segment, i.e. the stem), *strati* = layers, *raccolta* = collection, *confronta* = compare, *esegui* = run, *riproduci* = reproduce, *intestazione* = header.

## Main entry points

| file (`R/`) | main function | what it does |
|---|---|---|
| `sessione_pulita.R` | `carica_libreria_pulita()`, `riproduci_pulita()`, `confronta_manifesto()` | Loads the validated library from an explicit list of files in a clean R session, activates the faithful second-pass base (`define_base_forb_fedele`) and the Gauss-Newton cylinder solver (`least_squares_cylinder_fedele`), writes `manifesto_sessione_pulita.csv`; `riproduci_pulita()` reruns the reproduction of the earlier configuration; `confronta_manifesto()` compares the manifest with `docs/manifesto_validato.csv`. |
| `lancia_pulita.R` | (script) | Runs the clean-session reproduction from a terminal with `Rscript --vanilla`. |
| `riproduci_fedele.R` | `riproduci_fedele()` | Reruns `treeqsm` on three trees (3 seeds, both clouds) in a clean session and compares each quantity with the stored results of the final configuration. |
| `raccolta_modelli_quote.R` | `raccolta_modelli_quote()` | Runs `treeqsm` (8 parameter combinations per tree, several seeds) and collects, for each model, the stem diameter at the field heights 1.30, 2.00, 4.00 m (DBHcyl analogue) and the stem volumes. |
| `analisi_fedele.R` | `analisi_fedele()` | Summaries of the QSM diameters of the faithful port against the field data and against the stem-modelling procedure of Greco et al. (2023); paired tests and bootstrap. |
| `tabella7_fedele.R` | `tabella7_fedele()` | Tables 6 and 7 of the paper (errors at the three heights; paired comparisons direct fit, QSM, Greco et al.). |
| `confronta_octave_fedele.R` | `confronta_octave_fedele()` | Per-tree comparison of the R port with the original TreeQSM 2.4.1 run in GNU Octave (bias, limits of agreement, equivalence). |

## Direct estimator and extraction of the trees

| file | main function | what it does |
|---|---|---|
| `dbh_diretto.R` | `dbh_diretto()` (with `kasa_cerchio`, `hough_cerchio`, `cerchio_geometrico`) | Direct diameter at one height: 0.2 m slab, circle Hough transform, robust geometric least-squares fit. |
| `dbh_diretto_quote.R` | `valida_dbh_diretto_quote()` | The direct estimator at 1.30, 2.00 and 4.00 m, with recentring on the axis at each height. |
| `valida_dbh_diretto.R` | `valida_dbh_diretto()` | Validation of the direct estimator against the field data (seconds). |
| `estrazione_alberi.R` | `estrai_albero()` (with `carica_las_xyz`, `dtm_celle`, `rileva_tronchi`, `traccia_asse`) | Extraction of single trees from the plot LAS file: height above ground, stem detection, axis tracing, selection of the points. |
| `verifica_estrazione.R` | `verifica_estrazione()` | Checks the R extraction against the Python one. |

## Tests and diagnostics

| file | what it does |
|---|---|
| `esegui_test_r.R` (`esegui_test_r()`) | Runs the tests that do not need `treeqsm` (about one minute). |
| `prove_sintetiche.R` | Synthetic tests with generated stems, with the corrected code. |
| `diagnosi_errori_silenziosi.R` | Finds errors swallowed by `tryCatch` handlers during `treeqsm` (see Section 4 of the paper). |
| `diagnosi_stadi_R.R`, `profilo_fusto_R.R`, `strati_radice_R.R`, `fit_radice_R.R` | Intermediate quantities of the port at each stage, to be compared with the same quantities of the original in Octave (`octave/*.m`). |
| `prova_base_fedele.R`, `prova_gn_fedele.R`, `prova_tree_sets_fedele.R` | Tests of the faithful versions (orthonormal base of the second pass, Gauss-Newton solver, `tree_sets` functions) against the original in Octave. |
| `prova_orthonormal_casuale.R`, `variabilita_base_casuale.R` | Effect of the random orthonormal base used by the original on the results. |
| `correzioni_sessione_oggi.R` | Corrected versions of several functions of the port (name kept from the development session). |

## Other folders

| path | content |
|---|---|
| `octave/esegui_treeqsm_octave.m` | Runs the unmodified original TreeQSM in GNU Octave on the extracted trees. |
| `octave/diagnosi_stadi_octave.m`, `profilo_fusto_octave.m`, `strati_radice_octave.m`, `confronta_fit_radice_octave.m` | Counterparts in Octave of the diagnostic scripts in R. |
| `octave/GUIDA_OCTAVE.md` | Guide (in Italian) to the comparison with Octave. |
| `python/` | Auxiliary analyses: solver comparison, direct estimator, call graph. Paths must be adapted. |
| `docs/INTESTAZIONE_GPL.txt` | GPL notice added at the top of the R files (`add_gpl_headers.R`). |
| `docs/funzioni_inventario.csv`, `albero_chiamate.*` | Inventory of the functions and call tree of `treeqsm`. |
| `docs/manifesto_validato.csv`, `file_md5.txt`, `sessionInfo.txt` | Manifest of the validated library (function, defining file, md5 of the body), md5 of the files, session information. |

## Console messages

Messages printed by the scripts are in Italian. Control lines starting with `P_` (for example `P_C1`, `P_C5`) are checks written *before* the run and printed as TRUE/FALSE; `ATTENZIONE` means warning, `NON TROVATI` means not found.