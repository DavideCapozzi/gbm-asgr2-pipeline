# Research loop — hypothesis-centred (gbm-asgr2, v1, 2026-11-27)

The hypothesis drives every search. Each step must be able to change a decision; if the data or literature cannot support the hypothesis, the loop pivots instead of collecting more.

## Loop
**H** hypothesis + assumptions A1…An (each: tested / untested, cheapest test) → **L** literature, only for an untested assumption (flow §6 ladder, `scripts/lit/`) → **D** data, atlas/Census first, budget in `light_fetch.md` → **T** test with a deterministic script → **Gate** → continue · widen · pivot · stop → back to **H**.
State of the loop: `notes/LEDGER.md` (assumptions, cost spent, gate outcomes).

## Step contract (1–2 lines in the reply before the step; no extra tool call)
`Step <H|L|D|T> · tests A<n> · can change: <decision> · route: <cheapest> (alt: <easier one considered>) · cost: <calls, min, GB> · exit: <criterion>`
After the step: actual cost, and whether the decision moved. A step that cannot move a decision is skipped.

## Gates
- Criteria are numbers written **before** the data are seen (in the prompt or the LEDGER).
- **End-goal check** at every gate: does this still answer *which myeloid cells express which lectins, and can public data measure it?* (plus: is the lectin a usable target?)
- A pivot or a stop gets one blind `skeptic` pass (storyline only, no scores).

## Fallback ladder (cheapest first)
1. Same genes, better resolution: fine labels (Mo-/Mg-TAM, cDC2), deeper chemistry, protein (ADT).
2. Curated lectin list: lectins in general, ranked by prior evidence (`prior_counts.sh`) × detectability × subtype specificity × depth robustness × protein readout × targetability; Pace 2026 overlap is a bonus.
3. Module score of co-varying lectins: pools sparse genes into enough data points.
4. Reframe to protein: CITE-seq (c), GSM4972212.
5. Stop with a stated negative result (ADR): what was measured, why it is not enough, what data would be.

## Prompt feedback
Last line of every reply: `Flow: <evidence> → <change> · Prompt: <what in the prompt helped or misled> → <edit for next version>`. Nothing is written per reply; `/asgr-handoff` collects the lines into `notes/PROMPT_LOG.md` and applies ≤ 3 edits to the next prompt.

## Tools
- `scripts/lit/prior_counts.sh GENE…` — Europe PMC hit counts (glioma myeloid; + single-cell), one call for a list.
- `scripts/lit/facts.sh PMCID 'REGEX'` — streamed full-text snippets, nothing stored.
- `PILOT_GENES=<file> PILOT_ONLY=<id> scripts/pilot/census_detectability.py` — Census screen of any gene list (≈ 95 s, 10.8 GB per dataset).
- `scripts/pilot/pilot_checks.py --rank <cells.tsv.gz> --genes <file>` — `panel_rank.tsv` (detectability, donors, mac/mic, depth sign).
