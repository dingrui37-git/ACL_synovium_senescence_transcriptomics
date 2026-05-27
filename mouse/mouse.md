Mouse discovery minimal submission package
==========================================

Created: 2026-05-27 00:03:59 CEST
Package root: E:/R/ACLsenescence2/rebuild_submission/submission/mouse

This package was built by exact filename copying only.
No whole folder was copied. No fuzzy filename matching was used.
For each manifest entry, the organizer first checked the exact expected path.
If absent, it recursively searched only the two allowed roots using exact basename matching:
  1. E:/R/ACLsenescence2
  2. E:/R/ACLsenescence2 LD
When duplicate exact basenames were found, the preferred root and mouse-discovery/rebuild paths were prioritized, and all candidates were written to manifest/duplicate_exact_basename_candidates.csv.

Main contents:
  scripts/                       final ordered analysis scripts
  figures/                       manuscript Figure 1-3 image/PDF files
  tables_and_source_data/         exact source-data/result tables required by the figures and downstream derivations
  supplement/S1/                 supplemental sample-level MDS QC outputs
  logs_and_versions/             logs, sessionInfo, software-version records
  manifest/                      copy manifest, missing-file report, duplicate-candidate audit, package summary

Important files to inspect after running:
  manifest/package_summary.csv
  manifest/copy_manifest.csv
  manifest/MISSING_REQUIRED_FILES.csv, if present
  manifest/duplicate_exact_basename_candidates.csv, if present

If MISSING_REQUIRED_FILES.csv exists, the package is incomplete and should not be uploaded until those exact files are restored or their exact basenames are added to the script manifest.
