#!/usr/bin/env bash
# Light full-text check: stream a Europe PMC full text (or a local xml/txt), strip tags,
# print up to MAX snippets around REGEX. Stores nothing; HTTP error/timeout -> exit 2, no file written.
# Usage: scripts/lit/facts.sh PMC11328419 'CLEC12A|TGFBI' [context_chars=200]
#        MAX=10 scripts/lit/facts.sh literature/fulltext/x.xml 'ASGR2'
set -uo pipefail
export LC_ALL=C   # byte regex: .{0,200} under UTF-8 overflows some grep builds (ugrep)
src=$1; re=$2; n=${3:-200}
if [[ -f $src ]]; then txt=$(cat "$src")
else txt=$(curl -sf --max-time 30 "https://www.ebi.ac.uk/europepmc/webservices/rest/${src}/fullTextXML") \
       || { echo "FAIL $src (HTTP error or timeout; try OpenAlex pdf_url)" >&2; exit 2; }
fi
printf '%s' "$txt" | sed 's/<[^>]*>/ /g' | tr -s ' \n\t' ' ' | grep -oiE ".{0,$n}($re).{0,$n}" | head -n "${MAX:-5}"
exit 0
