#!/usr/bin/env bash
# Prior-evidence counts: Europe PMC hitCount per gene in title/abstract, glioma myeloid context.
# Usage: scripts/lit/prior_counts.sh ASGR2 SIGLEC9 ...   -> TSV: gene, glioma_myeloid, glioma_myeloid_sc
set -uo pipefail
U=https://www.ebi.ac.uk/europepmc/webservices/rest/search
hit() { curl -sf --max-time 30 -G "$U" --data-urlencode "query=$1" -d format=json -d resultType=idlist -d pageSize=1 \
        | grep -o '"hitCount":[0-9]*' | cut -d: -f2 || echo NA; }
printf 'gene\tglioma_myeloid\tglioma_myeloid_sc\n'
for g in "$@"; do
  q="TITLE_ABS:(glioma AND $g AND (macrophage OR microglia OR myeloid))"
  printf '%s\t%s\t%s\n' "$g" "$(hit "$q")" "$(hit "$q AND (\"single-cell\" OR scRNA-seq)")"
done
