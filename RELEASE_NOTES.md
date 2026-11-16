# Release v1.0.0 — code and result tables for the saphenous vein reproducibility manuscript

First public release of the analysis code, locked gene sets, archived result
tables, and reproduction evidence accompanying the manuscript
"Cross-cohort reproducibility of human saphenous vein remodeling signatures
and their directional alignment with occluded grafts".

## Contents

- `R/` — the full analysis chain (00 download + integrity verification,
  10-16 formal analyses and audits, 21-25 table/summary assembly, 17-20
  final figures), each script with hard semantic assertions.
- `gene_sets/` — the analysis plan locked on 2026-10-29, before the formal
  analyses were run.
- `results/` — the archived reference outputs (result tables, supplementary
  tables, figures).
- `reproduction_evidence/` — per-script logs, machine-readable audit
  (`reproduction_audit.csv`), 23 key-number/semantic checks
  (`key_number_checks.txt`), input inventory with SHA-256 for all 64 source
  inputs + 1 project manifest, byte-comparison notes, and session
  information from the verified clean-machine run.
- `run_reproduction_audit.R` — self-contained audit entry point: from the
  repository root it downloads all inputs from GEO, hash-locks them, runs
  the full chain, executes the checks, and compares every output
  bidirectionally against the archived reference; exits non-zero on any
  failure.
- `renv.lock` — locked package versions (R 4.4.2).

## How to reproduce

```r
Rscript run_reproduction_audit.R
```

Data files are downloaded from GEO at run time and are not redistributed.
See `README.md` for system prerequisites and the documented known limitation
(`renv::status()` metadata-only inconsistency for svglite/systemfonts/
textshaping on conda-based R installations).

## Verification status

All 16 scripts exit 0; 23/23 key checks PASS; 53/53 comparable outputs
byte-identical to the archived reference; 0 missing, 0 unexpected files.
