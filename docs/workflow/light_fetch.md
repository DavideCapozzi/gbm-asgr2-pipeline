# Light fetch — papers and datasets (gbm-asgr2, 2026-11-27)

General rules: `~/projects/flow.md` §6. This file = project weights, budget, measured tests.

## Weight and disk budget
w = importance (1 context · 2 supporting · 3 decides the storyline) × relevance (1 bulk/other tissue · 2 GBM sc · 3 GBM sc myeloid + our genes) × compatibility (0 unusable: mouse, restricted, odd format · 1 needs a new loader · 2 Census/existing loader) × curation (1 author files · 2 harmonised atlas).
Size is not in w: it sets the route. Big + low w → query remotely or drop.

| w | papers keep | datasets keep |
|---|---|---|
| 0 | nothing | nothing (note why) |
| 1–7 | id + 1-line fact in EVIDENCE | metadata ≤ 1 MB (sample sheet, rownames streamed) |
| 8–17 | + snippet in the fact (`scripts/lit/facts.sh`), no file | remote query (Census genes/obs) or processed subset ≤ 200 MB |
| ≥ 18 | + PDF/suppl ≤ 5 MB **only** if a figure/table value is cited | processed counts ≤ 2 GB, only the files the analysis reads; never FASTQ/BAM |

Project caps: full texts ≤ 20 MB, datasets ≤ 10 GB per machine. Before any download: `curl -sI` size and w; record id, w, bytes, decision fed.

## Current items
| item | w | kept | bytes |
|---|---|---|---|
| Pace 2026 suppl. (E1–E2, Table S1 values) | 3×3×2×1 = 18 | PDF | 1.5 MB |
| Dusoswa 2020 (E6–E9, figure values) | 3×2×1×1 = 6 → but figure values cited | PDF (exception, cited path) | 2.4 MB |
| 7 PMC full-text XMLs (E14–E24) | 4–12 | snippets/facts only; XMLs deleted (re-fetch with `facts.sh`) | 0 |
| GBmap Core/Extended (Census) | 3×3×2×2 = 36 | remote query, 14 genes | results only |
| GSE163120 human CSV + ADT | 3×3×1×1 = 9 (+ CITE for (c)) | processed CSVs ≈ 177 MB | 0 now |
| GSE182109 RAW.tar 2.46 GB | 2×2×1×1 = 4 if inside GBmap | drop if donors in GBmap; else processed only | 0 |
| GSE235676 / GSE154795 (1.1–1.3 GB) | 1×3×1×1 = 3 | metadata only | 0 |

## Tests (measured)
| test | old way | light way | result |
|---|---|---|---|
| T1 paper fact (E14, PMC11328419, CLEC12A/TGFBI BMD-TAM) | save fullTextXML (180,006 B) then grep | `scripts/lit/facts.sh PMC11328419 'CLEC12A'` streamed | same fact (snippet holds CLEC12A + TGFBI), 408 B kept (−99.8 %), 5 s (25 s when the network was shared with the pilot) |
| T1b failing full text (PMC11630528, HTTP 500) | 151 B error body saved as `.xml` (2 such files in `literature/fulltext/`) | `curl -sf` → exit 2, message, nothing written | no garbage files |
| T2 dataset, Core GBmap | h5ad download 8.13 GB (Extended 11.22 GB) | Census 14-gene query, one dataset per call: 86 s, 10.7 GB peak RAM, 0 B raw on disk, results 0.7 MB | −100 % raw disk; RAM is the cost → one dataset per query |
| T2b GSE182109 | RAW.tar 2.46 GB (`curl -sI`) | check donor overlap with GBmap from Census obs first: 0 B | download only if not in GBmap |

Gotcha found by the tests: interactive `grep` on the laptop is `ugrep`; `.{0,200}` under UTF-8 fails with "exceeds complexity limits" and returned 0 hits silently in a pipe. `facts.sh` runs non-interactive GNU grep with `LC_ALL=C`; a regex over-count check (Python `re`) gives the true number of mentions.

### Pilot (session 4)
- Glioma cohort atlas: 0 cells from Census 2025-11-08 in 25 s, 0.9 GB RAM → absence detected without a download.
- Only datasets that can change a decision were run (Core + replication); Extended/SL040/SL057/Darmanis skipped. The full 6-dataset run had stalled > 20 min on a too-small TileDB buffer.
