# Changelog

## 4.0.0 — 2026-08-17

- Added guarded IEDB BCR ingestion with paginated API downloads, explicit receptor/assay/epitope linkage, calculated-versus-curated chain annotations, and full-length nucleotide provenance.
- Added structured quantitative-measurement fields without unsafe unit normalization.
- Joined IEDB assay measurements through explicit receptor-group and assay IDs, with fail-closed support validation.
- Added lazy embedded PDB viewing with exact/homologous labels and RCSB fallback links.
- Added deterministic 99/95/90 heavy, light, and paired sequence-cluster indexes, family views, and redundancy exports.
- Added standalone and bundled BibTeX/RIS/CSV citation exports with source-aware deduplication.
- Added local paired `name,VH,VL` CSV/TSV screening, exact pair intersection, paired similarity metrics, and correct one-sequence-per-line parsing.
- Added conservative UniProt/HGNC/NCBI Taxonomy target identity metadata and explicit, validated protein/domain hierarchy with opt-in descendant search.
- Added dedicated entity summaries, stale-ID states, dated indexed-provenance timelines, and categorized correction reports.
- Added separate synthetic coding-DNA bundles with deterministic presets, AA digests, generation manifests, and strict ambiguous-residue rejection.
- Added an offline static-snapshot `pairs` CLI for target, antibody, exact-sequence, and indexed similarity queries with JSON/CSV/FASTA output.
- Added persistent cross-page selection, evidence-locked result lines, per-record FASTA, and selected-result ZIP bundles with metadata, citations, and reproducibility manifests.
- Added an explicit optional source-nucleotide contract and unavailable-state labelling; PAIRS does not back-translate amino-acid sequences implicitly.
- Added source-defensible CDR3/gene-call highlighting, selected-antibody comparison with global amino-acid alignments, and explicit source-record/exact-pair dataset modes.
- Added indexed full-chain VH/VL similarity search, aligned identity/coverage reporting, and related-sequence discovery from antibody records.
- Added conservative sequence-quality metadata, global source/origin/sequence/evidence filters, recent text views, and exportable local-only workspaces.
- Migrated the public static-data contract to schema version 4 under `data/v4/`.
- Separated direct-positive, functional-positive, and negative target page families.
- Replaced generic antibody target annotations with explicit direct, functional, negative, and literature collections.
- Added structure identity tiers, exact-structure-only facets, and PAIRS/upstream freshness metadata.
- Made fuzzy text matches suggestion-only; submitted searches navigate only for exact normalized resolutions.
- Scoped ordinary unpaired-chain identities to source records and kept confirmed VHH identity explicit.
- Quarantined malformed sequences; multispecific constructs now expose searchable sequence arms while unsupported arm-target claims remain quarantined.
- Added collision-safe unresolved target IDs, strain-preserving Pox targets, and labelled hierarchy derivations.
- Required complete source builds for production deployment and removed provenance-record truncation.
- Added verified Thera-SAbDab record deep links and explicit homepage-only provenance labels for sources without stable record routes.

## 3.0.0 — 2026-08-17

- Added local CDRH3/CDRL3 near-match search and batch FASTA screening.
- Added browse facets, V3 static-data access documentation, robots, and sitemap.
- Moved the public static-data contract to schema version 3 under `data/v3/`.

## 2.0.0 — 2026-08-16

- Renamed the project to **PAIRS — Pan-Antibody Integrated Retrieval System** across code, metadata, documentation, and download identifiers.
- Fixed PLAbDab semicolon-delimited literature co-mentions incorrectly becoming aliases of one target.
- Added regression tests proving unrelated co-mentioned targets remain separate.
- Added exact local sequence lookup for VH/VHH, VL, CDRH3, and CDRL3, including optional exact paired VH + VL matching.
- Added a browsable/sortable target catalogue.
- Added typo-tolerant suggestions and complete ArrowUp/ArrowDown/Enter keyboard navigation.
- Added direct `?ab=` antibody deep links.
- Added build-time pagination for large target result payloads.
- Added versioned `data/v2/` paths and a client-side schema compatibility guard.
- Added concurrent export page/shard loading.
- Added RCSB links for PDB identifiers.
- Added a visible partial-source warning banner.
- Added download retry/backoff and homepage discovery for every enabled V2 source.
- Added cross-shard referential-integrity validation.
- Added human-review target-alias candidate output without automatic merging.
- Added dark-mode variables and retained horizontally scrollable example chips on small screens.
- Added OpenGraph/Twitter metadata and a social preview image.
- Added Ruff and Prettier CI gates.

## 1.0.0 — 2026-08-16

- Initial static GitHub Pages application with no runtime backend.
- Public-source adapters for PLAbDab, Thera-SAbDab, CoV-AbDab and Pox-AbDab.
- Exact VH/VL pair deduplication and provenance-preserving interaction model.
