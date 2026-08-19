# Architecture

```text
public source homepages
        ↓ discover current download link
source downloads + retry/backoff
        ↓
source adapters
        ↓
provenance-preserving observations
        ↓
conservative normalization
        ↓
exact sequence-pair deduplication
        ↓
target + interaction model
        ↓
versioned static indexes / paged shards
        ↓
GitHub Pages / any static HTTP server
```

There is no runtime backend. Python exists only at build time. V4 writes exact-chain lookup files and deduplicated CDRH3/CDRL3 length buckets; browser near matching fetches only the query length and adjacent buckets before applying bounded Levenshtein. Approximate VH/VL matching is intentionally deferred because it requires alignment-aware retrieval. Target retrieval is semantic: direct-positive observations are the default, while functional and negative observations have explicit page families and literature context never creates target hits.

## Static data layout

- `data/v4/manifest.json` — schema/app version, PAIRS snapshot date, upstream freshness, source status, and counts
- `data/v4/targets.json` — compact searchable/browsable target index
- `data/v4/targets/<hash>/page-NNN.json` — direct-positive paginated results for one target
- `data/v4/targets/<hash>/functional-page-NNN.json` — functional-positive observations for one target
- `data/v4/targets/<hash>/negative-page-NNN.json` — negative observations for one target
- `data/v4/targets/<hash>/index.json` — page manifest for all three target result families
- `data/v4/antibody-search/xx.json` — compact two-character name/alias buckets for autocomplete and fuzzy suggestions
- `data/v4/antibody-exact/xx.json` — lossless two-character buckets mapping normalized antibody names/aliases directly to antibody IDs
- `data/v4/sequence-search/xx.json` — SHA-256-prefixed exact sequence lookup buckets
- `data/v4/antibodies/xx.json` — full antibody records keyed by content-derived antibody ID
- `data/v4/target-review.json` — human-review alias candidates from direct-positive target assignments; never consumed as automatic merges

## Runtime

The homepage first reads the V4 manifest and refuses to proceed when the schema is incompatible. It then loads only `targets.json`. Selecting a target fetches direct-positive page one; functional and negative pages load only after the corresponding explicit filter is selected. Changing to a global non-default sort intentionally loads all pages for the active evidence family.

Antibody-name search fetches only the relevant two-character prefix bucket. Exact textual lookup uses the dedicated lossless `antibody-exact` map; fuzzy suggestions use the compact `antibody-search` records and never resolve a submitted query automatically. Text terms are Unicode NFKD-normalized, case-folded, converted to ASCII alphanumeric runs separated by single spaces, and bucketed from their compact first two characters. This intentionally treats punctuation variants such as `PD-L1` and `PDL1` as equivalent while leaving sequence normalization separate. Exact lookup may return multiple antibody IDs, which the UI presents for user selection. Exact sequence search hashes the normalized amino-acid sequence in the browser and fetches only the corresponding two-character SHA-256 bucket. Full antibody records stay in deterministic two-character shards and load only for detail views/exports. Their semantic collections are explicit: `direct_targets`, `functional_targets`, `negative_evidence`, and `literature_mentions`. Multispecific sequence arms retain construct membership without inheriting construct-level targets. Structure records retain `structure_tiers`, and UI structure facets count exact (100% identity) records only. The manifest distinguishes the PAIRS indexed date from each source's upstream `last_modified` value; provenance distinguishes verified record URLs from source-homepage links.

No search term or sequence needs to leave the browser.
