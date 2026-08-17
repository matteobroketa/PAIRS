# PAIRS static data access

PAIRS has no runtime API server. Immutable static JSON files are published for each snapshot; this is a read-only data interface, not an SLA-backed REST API.

V3 uses `schema_version: 3`. Start with `data/v3/manifest.json` and `data/v3/targets.json`. Target result pages are at `data/v3/targets/<dir>/page-NNN.json`, antibody records use `data/v3/antibodies/<id-prefix>.json`, exact lookup uses SHA-256 prefix files in `data/v3/sequence-search/`, and CDRH3/CDRL3 candidate buckets are at `data/v3/sequence/cdrh3/NN.json` and `cdrl3/NN.json`.

```sh
curl "$PAIRS_BASE/data/v3/manifest.json"
```

```python
import requests
targets = requests.get(f"{base}/data/v3/targets.json").json()
```

```js
const manifest = await fetch(`${base}/data/v3/manifest.json`).then(r => r.json());
```
