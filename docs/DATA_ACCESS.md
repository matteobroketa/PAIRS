# PAIRS static data access

PAIRS has no runtime API server. Immutable static JSON files are published for each snapshot; this is a read-only data interface, not an SLA-backed REST API.

V4 uses `schema_version: 4`. Start with `data/v4/manifest.json` and `data/v4/targets.json`. Direct-positive target result pages are at `data/v4/targets/<dir>/page-NNN.json`; functional-positive and negative observations are separate `functional-page-NNN.json` and `negative-page-NNN.json` families. Antibody records use `data/v4/antibodies/<id-prefix>.json`, exact lookup uses SHA-256 prefix files in `data/v4/sequence-search/`, and CDRH3/CDRL3 candidate buckets are at `data/v4/sequence/cdrh3/NN.json` and `cdrl3/NN.json`.

Full-chain candidate indexes are under `data/v4/similarity/heavy/` and `light/`. Their `index.json` files version the retrieval contract. PAIRS uses bottom-k hashes of distinct amino-acid 5-mers only to retrieve candidates, then ranks displayed results using deterministic global amino-acid alignment. A signature hit is not itself a similarity result, and sequence similarity is not evidence of shared specificity or lineage.

The default target page represents only direct-positive `binds`/`targets` relationships. Functional relationships (`neutralizes`/`protects`) and negative relationships are opt-in page families. Literature context (`mentioned_with`) is retained in antibody `literature_mentions` collections but never creates a target result page. Antibody records use explicit `direct_targets`, `functional_targets`, `negative_evidence`, and `literature_mentions` collections; consumers must not infer semantics from a generic `targets` field. Structure records are tiered through `structure_tiers`; only the 100% identity tier is an exact structure.

Multispecific Thera-SAbDab records expose `constructs` and `arms` so each available sequence arm remains searchable. Construct-level targets never become arm targets unless an explicit upstream arm-target field supports the assignment. Provenance uses `record_url` with `link_scope: "record"` only for verified deep links; `link_scope: "source_homepage"` is not record-level verification.

Optional `vh_nt_source` and `vl_nt_source` fields contain only source-reported nucleotide sequences and require matching `nucleotide_provenance`. Missing fields mean unavailable source DNA; consumers must not infer nucleotide sequence from amino-acid fields. Browser export bundles contain machine-readable AA FASTA, separate heavy/light FASTA, CDR3 FASTA, metadata CSV, citations, and a reproducibility manifest. A source-nucleotide FASTA is included only when source data exist.

IEDB-derived antibody records may expose `chain_annotations` with separate `calculated` and `curated` region/gene calls. Interaction records may include `target_external_id`, `assay_ids`, `receptor_group_id`, and structured `measurements`; consumers must preserve the original unit and qualifier.

Dataset exports state their redundancy mode explicitly. `source_records` repeats selected sequence entities per upstream record, while `unique_exact_vh_vl_pairs` deduplicates only complete normalized VH+VL pairs and keeps incomplete chains source-scoped. The 99/95/90 modes use indexed candidate retrieval followed by global edit identity, at least 90% length coverage, and single-link clustering; paired records must pass VH and VL separately. Candidate retrieval can have false negatives. These are sequence clusters, not clonal lineages or specificity claims.

Structure buttons embed Mol* only after a user click and only for syntactically valid four-character PDB IDs. Exact and homologous tiers remain distinct. Structure metadata is sequence-level and does not assert an antibody–target complex; every view retains an RCSB fallback link.

Selected-result bundles include deduplicated `references.bib`, `references.ris`, and `references.csv` files derived only from available upstream reference strings, source record IDs, and links. No missing bibliographic fields are invented.

Batch screening accepts multi-FASTA, one amino-acid sequence per line, or RFC-4180-style CSV/TSV with `name,VH,VL` columns. Paired matches require the VH and VL to resolve to the same PAIRS sequence entity; rows are never paired by result order. Similarity output reports VH and VL identity/coverage separately.

## Local CLI

Install the repository locally with `python -m pip install -e .`, then query a downloaded static snapshot without a server or network access:

```sh
pairs --data data/v4 --format fasta target HER2 --paired
pairs --data data/v4 target "SARS-CoV-2 Spike" --include-descendants
pairs --data data/v4 sequence EVQLVESGG --chain heavy
```

Target searches are direct-positive and exact-scope by default. Functional, negative, and descendant records require explicit flags. JSON output includes snapshot, query semantics, and PAIRS IDs.

## Generated coding DNA

`Generate coding DNA bundle` is an explicit selected-record operation. It creates `GENERATED_FROM_AA.fasta`, a generation manifest with preset/version and input AA SHA-256, and a warning README. Generated DNA is a deterministic variable-domain back-translation—not source DNA, expression optimization, or a complete construct—and is never written into `vh_nt_source`, `vl_nt_source`, or source-nucleotide exports.

Browser workspaces are local UI state, not part of the published data API. Workspace JSON contains selected PAIRS IDs, recent target/antibody views, current filter settings, and snapshot identity. Raw sequence queries are never stored. Imports retain IDs from older snapshots but show a mismatch warning rather than silently reinterpreting them.

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
