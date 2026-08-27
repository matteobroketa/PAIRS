# Data sources

PAIRS is an aggregation/search layer. It does **not** claim ownership of imported scientific data. Every result preserves source provenance, and researchers should cite the originating database/publication.

## Enabled in the V4 public build

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

### IEDB BCR

The Immune Epitope Database BCR export is ingested through its documented query API under CC BY 4.0. PAIRS preserves any explicit receptor-group/assay/epitope identifiers and keeps calculated and curated chain annotations separate. Calculated variable domains are searchable; source nucleotide strings remain labelled as full-length BCR sequences. The current export is sequence-oriented and its epitope columns may be empty, so rows without explicit target linkage create no interaction.

Upstream: https://www.iedb.org/database_export_v3.php

### SAbDab2

The Structural Antibody Database (OPIG) is ingested from its public API and
bulk instance summary under CC BY 4.0. Canonical IMGT-numbered variable-domain
sequences are joined to explicit antibody-instance/antigen pairings. PDB
structures are retained as provenance, and only those explicit pairings create
`STRUCTURE`/`binds` observations; free-text titles and keywords are not target
evidence. The downloader stores the API pages and summary together in the
local source cache so the public build does not need to redistribute upstream
structure files.

Upstream: https://sabdab.opig.stats.ox.ac.uk/

## Cross-source normalization and reconciliation

Target resolution prefers source-provided UniProt, HGNC, and NCBI Gene identifiers. The
curated `config/target_aliases.json` file is used for reviewed synonyms and overrides;
unresolved stable identifiers remain distinct rather than being collapsed by a free-text slug.
SAbDab2 antigen identifiers and antigen organism context are retained on each interaction.

Antibody observations expose an explicit sequence architecture (`paired_vh_vl`, `heavy_only`,
`light_only`, `vhh`, `scfv`, `cdr_only`, `nucleotide_only`, or `sequence_unavailable`). Records
are reconciled across sources only with strong evidence such as identical variable sequences,
matching stable sequence accessions, explicit source aliases, or shared structural provenance.
Each contributing source record and interaction link remains attached after reconciliation.

When a curated source supplies a VH/VL GenBank, ENA, or EMBL accession, the build retrieves only
that exact record, caches the flatfile under `.cache/sources/accessions/`, and stores its URL,
accession, checksum, and chain scope. PAIRS does not crawl GenBank or ENA globally.

## Download freshness

The four legacy OPIG CSV sources use homepage discovery in `config/sources.json`. SAbDab2 uses deterministic API pagination plus its all-summary CSV; IEDB uses deterministic PostgREST CSV pagination against `bcr_export` plus exact receptor-group/assay joins from `bcr_search` to `bcell_export`. No epitope-name join is used. Download or support-schema failures stop complete builds; downloads use bounded retry/backoff.

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
