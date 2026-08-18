# Biological gold standard

`tests/gold/cases.json` is the permanent registry of manually reviewed real source records. Each case
pins the upstream source record, review date, PAIRS sequence entity, chain state, relationship family,
target interpretation, structure tier, and provenance expectation. `tests/test_gold_records.py`
checks those claims against the versioned V4 starter snapshot.

The initial registry is deliberately small and covers Thera-SAbDab, CoV-AbDab, Pox-AbDab, and
PLAbDab. Existing IEDB fixture tests independently check receptor-group/assay measurement linkage.
This seed is infrastructure, not sufficient empirical validation. Expansion to 200–500 reviewed
records—especially natural relatives, aliases, negative assays, hierarchies, multispecifics, and IEDB
measurements—remains a curation task. Cases must be added only after checking the cited source record;
generated or inferred examples must not be labeled gold.
