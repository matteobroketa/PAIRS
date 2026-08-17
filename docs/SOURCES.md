# Data sources

PAIRS is an aggregation/search layer. It does **not** claim ownership of imported scientific data. Every result preserves source provenance, and researchers should cite the originating database/publication.

## Enabled in the V2 public build

### PLAbDab

Patent & Literature Antibody Database (Oxford Protein Informatics Group). The importer uses the public paired-sequence download and retains `targets_mentioned` as **literature metadata**, not direct binding evidence. Semicolon-delimited co-mentioned target terms are emitted as independent `mentioned_with` interactions and are never made aliases of one another.

Upstream: https://opig.stats.ox.ac.uk/webapps/plabdab/

### Thera-SAbDab

Therapeutic Structural Antibody Database (OPIG). Used for therapeutic sequence, target, format, status, and structural metadata.

Upstream: https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/therasabdab/search/

### CoV-AbDab

Coronavirus-binding antibody database (OPIG), distributed upstream under CC BY 4.0. PAIRS preserves positive and negative binding/neutralisation annotations.

Upstream: https://opig.stats.ox.ac.uk/webapps/covabdab/

### Pox-AbDab

Orthopoxvirus antibody database (OPIG), distributed upstream under CC BY 4.0. PAIRS preserves binding, neutralisation, and protection annotations, including negative evidence.

Upstream: https://opig.stats.ox.ac.uk/webapps/poxabdab/

## Download freshness

All four enabled sources use homepage discovery in `config/sources.json`. If an upstream file name rotates and the configured discovery expression no longer finds the current download, the build marks that source failed rather than silently continuing to rely on an old hard-coded file path. Downloads use bounded retry/backoff.

## Deliberately disabled by default

### ABCD

ABCD is highly relevant, but its documented full-database route involves contacting the maintainers. The repository reserves an integration slot but does not scrape or redistribute ABCD by default.

### NaturalAntibody ASD / Patents

These are scientifically valuable sources whose publicly described data-access conditions include non-commercial research restrictions. They are not embedded in the default public GitHub Pages build. A future local/private adapter can be added without changing the core schema.

## Adding a source

1. Add a source definition to `config/sources.json`.
2. Implement a generator in `pipeline/sources.py` returning an `AntibodyObservation` plus one or more `InteractionObservation` objects.
3. Register it in `ADAPTERS`.
4. Add fixture tests for its column mapping and relationship semantics.
5. Document license, attribution, download discovery, and evidence meaning here.

Never silently convert a literature keyword into `binds`. Use the least-strong relationship/evidence class justified by the source.
