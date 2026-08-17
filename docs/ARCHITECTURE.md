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

There is no runtime backend. Python exists only at build time. V3 writes exact-chain lookup files and deduplicated CDRH3/CDRL3 length buckets; browser near matching fetches only the query length and adjacent buckets before applying bounded Levenshtein. Approximate VH/VL matching is intentionally deferred because it requires alignment-aware retrieval.

## Static data layout

- `data/v2/manifest.json` — schema/app version, snapshot, source status, and counts
- `data/v2/targets.json` — compact searchable/browsable target index
- `data/v2/targets/<hash>/page-NNN.json` — build-time paginated results for one target
- `data/v2/targets/<hash>/index.json` — page manifest for a target
- `data/v2/antibody-search/xx.json` — two-character name/alias prefix buckets
- `data/v2/sequence-search/xx.json` — SHA-256-prefixed exact sequence lookup buckets
- `data/v2/antibodies/xx.json` — full antibody records keyed by content-derived antibody ID
- `data/v2/target-review.json` — human-review alias candidates; never consumed as automatic merges

## Runtime

The homepage first reads the V2 manifest and refuses to proceed when the schema is incompatible. It then loads only `targets.json`. Selecting a target fetches page one; additional target pages load on demand. Changing to a global non-default sort intentionally loads all pages for that target.

Antibody-name search fetches only the relevant two-character prefix bucket. Exact sequence search hashes the normalized amino-acid sequence in the browser and fetches only the corresponding two-character SHA-256 bucket. Full antibody records stay in deterministic two-character shards and load only for detail views/exports.

No search term or sequence needs to leave the browser.
