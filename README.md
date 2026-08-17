# PAIRS — Pan-Antibody Integrated Retrieval System

PAIRS is a static, provenance-aware search interface for public antibody sequences, targets, structures, and source-level evidence.

**V2 public sources:** PLAbDab, Thera-SAbDab, CoV-AbDab, and Pox-AbDab. The public build deliberately excludes sources whose bulk redistribution terms require additional permission or non-commercial gating.

## Deploy to GitHub Pages

1. Create an empty repository and copy this entire folder into it.
2. Commit and push to `main`.
3. In **Settings → Pages**, set **Source** to **GitHub Actions**.
4. The included `deploy-pages.yml` discovers the current download links for enabled sources, downloads them with retry/backoff, builds the full V2 static index, validates cross-file integrity, and publishes the site.

No runtime server, database, API key, Node build step, or paid service is required. Node is used only in CI/development to enforce Prettier formatting.

## Run locally

The tracked `data/v2/` directory contains a compact starter snapshot so the UI works immediately over HTTP:

```bash
python -m http.server 8000
```

Then open `http://localhost:8000`.

Do not open `index.html` via `file://`; browsers block the JSON fetches used by the app.

## Build the current public index locally

```bash
python -m pipeline.build --output data/v2
python -m pipeline.validate data/v2
```

Python 3.10+ is sufficient; the data pipeline uses the standard library only.

The tracked starter snapshot uses a 750-record-per-source cap (sources with fewer records contribute all available rows). A fresh comparable demo build is:

```bash
python -m pipeline.build --output data/v2 --max-records 750
```

To reuse already-downloaded source files:

```bash
python -m pipeline.build --output data/v2 --offline --cache .cache/sources --max-records 750
```

## What V2 does

- one search surface for normalized targets and antibody/therapeutic names
- exact browser-side VH/VHH, VL, CDRH3, and CDRL3 sequence lookup
- optional exact paired VH + VL lookup
- browsable and sortable target catalogue
- typo-tolerant target and antibody suggestions
- keyboard-accessible search suggestions
- exact sequence-pair deduplication
- separate target, antibody, and interaction entities
- positive **and negative** binding/neutralisation/protection evidence
- `CURATED` and `LITERATURE_METADATA` kept distinct; no opaque confidence score
- PLAbDab multi-target literature co-mentions fanned into independent `mentioned_with` interactions
- build-time pagination of large target payloads
- lazy browser-only target pages and antibody shards
- shareable target URLs and direct `?ab=` antibody links
- RCSB links for PDB structures
- CSV and FASTA export with concurrent shard loading
- source-status banner when a scheduled build is partial
- versioned `data/v2/` path plus client-side schema guard
- build-time alias-review candidates based on overlapping antibody sets; candidates are never auto-merged
- weekly automated rebuild workflow with download discovery and retry/backoff
- no runtime calls to third-party antibody databases
- no analytics or sequence uploads

## Target correctness policy

PAIRS does not infer synonyms merely because target names occur in the same paper. In particular, PLAbDab's semicolon-delimited `targets_mentioned` field is split into independent target interactions before normalization. Co-occurrence never adds one target name to another target's alias set.

`config/target_aliases.json` remains deliberately conservative. V2 also generates `data/v2/target-review.json`, a human-review queue of target pairs with heavily overlapping antibody sets. This is advisory curation data only; the build never merges those candidates automatically.

The next normalization step should be authoritative identifier mapping at build time (for example HGNC/UniProt for human protein targets and NCBI Taxonomy for pathogens) while retaining unresolved text targets rather than guessing.

## Scientific caution

PAIRS is a search/indexing tool, not an experimental validation authority. PLAbDab `targets_mentioned` records are represented as `LITERATURE_METADATA`, not measured binding. Source databases can contain incomplete sequences, naming ambiguity, legacy annotations, or publication-derived metadata. Verify important records against the originating source before experimental use.

## Data and licensing

Repository software: MIT. Imported scientific data: upstream terms apply. See `docs/SOURCES.md` and `docs/DATA_LICENSING.md`.

## Development

```bash
python -m pip install pytest ruff
pytest -q
ruff check pipeline tests
npx --yes prettier@3.9.6 --check index.html assets config
python -m pipeline.validate data/v2
```

`make test`, `make lint`, `make demo`, and `make serve` wrap the common commands.

## Repository map

```text
index.html                 static app shell
assets/                    browser JS/CSS + social preview image
pipeline/                  Python import/index compiler
config/                    sources + conservative target aliases
data/v2/                   compact tracked starter snapshot
.github/workflows/         tests, Pages deployment, scheduled rebuild
docs/                      architecture, schema, sources and licensing
```

## Citation

If this project becomes part of published work, cite PAIRS **and** the source databases contributing the records used in the analysis. Upstream citations are listed in `config/sources.json` and `docs/SOURCES.md`.
