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

### NanoSaurus

NanoSaurus provides public single-domain antibody records with explicit antigen
links and UniProt/PDB identifiers. PAIRS retrieves its public JSON query pages
deterministically, classifies records as VHH, and preserves the entry URL and
validation metadata. The build does not infer targets from nanobody names.

Upstream: https://www.nanosaurus.org/

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
Recovered nucleotide fields are also written to the separate `nucleotide-search` exact index;
they are never translated implicitly into amino-acid sequences.

## Download freshness

The four legacy OPIG CSV sources use homepage discovery in `config/sources.json`. SAbDab2 uses deterministic API pagination plus its all-summary CSV; IEDB uses deterministic PostgREST CSV pagination against `bcr_export` plus exact receptor-group/assay joins from `bcr_search` to `bcell_export`. No epitope-name join is used. Download or support-schema failures stop complete builds; downloads use bounded retry/backoff.

## Deliberately disabled by default

### ABCD

ABCD is highly relevant, but its documented full-database route involves contacting the maintainers. The repository reserves an integration slot but does not scrape or redistribute ABCD by default.

### NaturalAntibody ASD / Patents

These are scientifically valuable sources whose publicly described data-access conditions include non-commercial research restrictions. They are not embedded in the default public GitHub Pages build.

PAIRS includes an offline ASD adapter for an authorized CSV or JSON export. It
accepts the published `heavy_sequence`, `light_sequence`, `affinity_type`,
`affinity`, `antigen_sequence`, `metadata.target_uniprot`,
`metadata.target_pdb`, `metadata.target_name`, `metadata.source_url`,
`confidence`, and `scfv` fields. Stable antigen identifiers are preferred over
free-text names, and measurements/source URLs remain attached to each ASD row.
The adapter is intentionally not downloadable and is not selected by the public
build. A local user can run `python -m pipeline.build --source
naturalantibody_asd --offline --cache PATH --output PATH` against a legitimately
obtained export. Do not commit or publish the restricted export or its generated
data.

### IMGT/mAb-DB

IMGT/mAb-DB is a curated therapeutic-antibody resource with specificity and
format metadata. Its documented interface supports query/result exports, but
PAIRS does not crawl the query site. An authorized CSV/JSON export can be
loaded with the disabled-by-default `imgt_mabdb` adapter; target identifiers,
format, status, structures, and record links are preserved.

## Remaining coverage gaps

The follow-up source audit found several relevant resources that are not direct
public-build inputs yet:

- OAS is a large, openly licensed antibody-sequence repository, but its public
  metadata is generally repertoire/study context rather than antibody-target
  interaction evidence. The optional `oas` adapter accepts explicitly selected
  data-unit CSV/CSV.GZ files as sequence-only records and creates no target
  edges without an explicit target-linkage feed.
- IMGT/LIGM-DB is valuable for sequence recovery and immunoglobulin annotation,
  but it is not a general antibody-antigen interaction database. It remains a
  future sequence-enrichment source, not a target-claim source.
- AntibodyDB exposes a relevant influenza collection, but its currently
  documented route requires account/request approval. It should be integrated
  only after an authorized export/API route and terms for redistribution are
  confirmed.
- ABCD's documented full-database route also requires maintainer contact, so it
  remains an offline authorized-export adapter rather than a scraper.

These sources should not be approximated by extracting names from search pages:
PAIRS needs a stable target identifier or an explicitly curated source
relationship before creating a target edge.

## Adding a source

1. Add a source definition to `config/sources.json`.
2. Implement a generator in `pipeline/sources.py` returning an `AntibodyObservation` plus one or more `InteractionObservation` objects.
3. Register it in `ADAPTERS`.
4. Add fixture tests for its column mapping and relationship semantics.
5. Document license, attribution, download discovery, and evidence meaning here.

Never silently convert a literature keyword into `binds`. Use the least-strong relationship/evidence class justified by the source.
