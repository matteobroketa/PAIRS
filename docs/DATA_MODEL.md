# Data model

PAIRS deliberately separates **antibodies**, **targets**, and **source-level interactions**.

## Antibody

An antibody is automatically deduplicated only when its normalized heavy/light sequence pair is identical. Heavy-only, light-only, and unsequenced records remain distinct unless an exact sequence identity rule applies. Names alone never trigger a merge.

Core fields include stable content-derived ID, preferred name, aliases, VH/VHH, VL, CDRH3/CDRL3 when supplied, V/J calls, format, organism/origin, structures, therapeutic status, and source records.

## Target

Targets are conservatively normalized using `config/target_aliases.json`. V2 uses deterministic local IDs (`target:<slug>`) rather than pretending every text term has been confidently mapped to an authoritative biological identifier.

Co-occurrence is **not** synonymy. A semicolon-delimited PLAbDab literature mention group is fanned into separate interactions before target resolution. The build never inserts one co-mentioned target term into another target's alias set.

`data/v2/target-review.json` reports target pairs whose antibody sets overlap strongly so a human can inspect possible missed aliases. Those suggestions never change the live index automatically.

A future identifier layer should map resolvable protein targets to authoritative IDs (for example HGNC/UniProt) and pathogens to NCBI Taxonomy while preserving unresolved text labels rather than guessing.

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

This prevents destructive flattening of positive/negative or experimental/metadata evidence.

## Evidence classes

V2 currently emits `CURATED` and `LITERATURE_METADATA`. The schema reserves stronger/future categories such as `STRUCTURE`, `MEASURED`, `PUBLICATION`, `PATENT`, `METADATA`, and `SIMILARITY`. The UI intentionally does not present an opaque numeric confidence score.
