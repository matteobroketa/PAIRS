# Data model

PAIRS deliberately separates **antibodies**, **targets**, and **source-level interactions**.

## Antibody

An antibody is automatically deduplicated only when its normalized heavy/light sequence pair is identical. Heavy-only, light-only, and unsequenced records remain distinct unless an exact sequence identity rule applies. Names alone never trigger a merge.

Core fields include stable content-derived ID, preferred name, aliases, VH/VHH, VL, source-provided CDR1/2/3 regions, V/D/J calls, format, organism/origin, structures, therapeutic status, and source records. `chain_annotations` retains calculated and curated calls separately with their source coordinates and annotation origin.

Source-reported nucleotide sequences use the optional `vh_nt_source` and `vl_nt_source` fields with chain-specific `nucleotide_provenance`. Their absence means no source nucleotide sequence is available. PAIRS never treats a back-translation as source DNA.

`sequence_quality` contains only derivable or source-reported facts: chain pairing/availability, chain lengths, explicitly labelled VHH status, ambiguous/non-standard residues, quarantine state, native nucleotide availability, and source-reported completeness. Length alone never makes a sequence “complete”; absent source evidence is represented as `unknown_not_inferred`.

## Target

Targets are conservatively normalized using `config/target_aliases.json`. V4 uses deterministic local IDs (`target:<slug>`) rather than pretending every text term has been confidently mapped to an authoritative biological identifier.

Co-occurrence is **not** synonymy. A semicolon-delimited PLAbDab literature mention group is fanned into separate interactions before target resolution. The build never inserts one co-mentioned target term into another target's alias set.

`data/v4/target-review.json` reports target pairs whose antibody sets overlap strongly so a human can inspect possible missed aliases. Those suggestions never change the live index automatically.

The identifier layer maps only resolvable protein targets to authoritative IDs (for example HGNC/UniProt) and selected pathogens to NCBI Taxonomy while preserving unresolved text labels rather than guessing.

V4 now provides a conservative first identifier layer in `entity`: verified identifiers supplement, but never replace, stable PAIRS target IDs. `hierarchy.parent_id`, `relation`, and reciprocal `children` express explicitly curated protein/domain containment. Every mapping records its scope and verification provenance; hierarchy expansion is never implicit.

## Interaction

Every source assertion becomes an interaction with:

- antibody ID
- target ID and original target text
- relationship (`binds`, `does_not_bind`, `neutralizes`, `does_not_neutralize`, `protects`, `targets`, `mentioned_with`, …)
- evidence class
- source + source record ID
- reference
- epitope/domain annotation
- note/qualifier

Optional structured `measurements` retain metric, numeric value when parseable, qualifier, original value, unit, and assay context without cross-unit normalization. `target_external_id`, `assay_ids`, and `receptor_group_id` preserve explicit upstream linkage.

This prevents destructive flattening of positive/negative or experimental/metadata evidence.

Source records may carry `record_date` and `record_date_field` only when the upstream record supplies an applicable date. The UI orders these as an “indexed provenance timeline”; it is not a claim about invention, priority, or first publication.

## Evidence classes

V4 emits `CURATED`, `MEASURED`, and `LITERATURE_METADATA` where the source contract supports them. The schema reserves additional categories such as `STRUCTURE`, `PUBLICATION`, `PATENT`, `METADATA`, and `SIMILARITY`. The UI intentionally does not present an opaque numeric confidence score.
