# Target identity and hierarchy

PAIRS keeps stable local target IDs even when an authoritative external accession is available. `config/target_entities.json` adds manually verified entity type, organism, UniProt/HGNC/NCBI Taxonomy identifiers, and explicit parent relations to a conservative subset of unambiguous canonical targets.

External identifiers never cause fuzzy merges. Unresolved source terms retain collision-safe local IDs, and domains remain distinct targets. Hierarchy edges currently use `protein_of` or `subdomain_of`, include mapping provenance, and are validated for missing parents, non-reciprocal child edges, and cycles.

Target searches default to exact-target results. `Include subdomains/children` is an explicit opt-in and does not rewrite child evidence as direct evidence for the parent.
