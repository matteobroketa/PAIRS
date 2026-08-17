# PAIRS static data access

PAIRS has no runtime API server. Immutable static JSON files are published for each snapshot; this is a read-only data interface, not an SLA-backed REST API.

V4 uses `schema_version: 4`. Start with `data/v4/manifest.json` and `data/v4/targets.json`. Direct-positive target result pages are at `data/v4/targets/<dir>/page-NNN.json`; functional-positive and negative observations are separate `functional-page-NNN.json` and `negative-page-NNN.json` families. Antibody records use `data/v4/antibodies/<id-prefix>.json`, exact lookup uses SHA-256 prefix files in `data/v4/sequence-search/`, and CDRH3/CDRL3 candidate buckets are at `data/v4/sequence/cdrh3/NN.json` and `cdrl3/NN.json`.

The default target page represents only direct-positive `binds`/`targets` relationships. Functional relationships (`neutralizes`/`protects`) and negative relationships are opt-in page families. Literature context (`mentioned_with`) is retained in antibody `literature_mentions` collections but never creates a target result page. Antibody records use explicit `direct_targets`, `functional_targets`, `negative_evidence`, and `literature_mentions` collections; consumers must not infer semantics from a generic `targets` field. Structure records are tiered through `structure_tiers`; only the 100% identity tier is an exact structure.

Multispecific Thera-SAbDab records expose `constructs` and `arms` so each available sequence arm remains searchable. Construct-level targets never become arm targets unless an explicit upstream arm-target field supports the assignment. Provenance uses `record_url` with `link_scope: "record"` only for verified deep links; `link_scope: "source_homepage"` is not record-level verification.

```sh
curl "$PAIRS_BASE/data/v4/manifest.json"
```

```python
import requests
targets = requests.get(f"{base}/data/v4/targets.json").json()
```

```js
const manifest = await fetch(`${base}/data/v4/manifest.json`).then(r => r.json());
```
