# Contributing

Contributions are welcome, especially new public data-source adapters, authoritative target identifiers, target alias corrections, provenance improvements, and reproducible bug reports.

## Development

```bash
python -m pip install pytest ruff
pytest -q
ruff check pipeline tests
npx --yes prettier@3.9.6 --check index.html assets config
python -m pipeline.build --output data/v4 --max-records 750
python -m pipeline.validate data/v4
python -m http.server 8000
```

The data pipeline intentionally has no third-party runtime dependencies. Ruff and Prettier are development/CI checks only.

## Source adapters

A source contribution must include:

1. Current upstream homepage and machine-readable access route.
2. License/redistribution note and preferred citation.
3. A source adapter returning `AntibodyObservation` and `InteractionObservation` records.
4. Required-column/schema checks.
5. Tests demonstrating at least one positive record and, when available, negative evidence handling.
6. Documentation in `docs/SOURCES.md`.

Do not upgrade evidence semantics. A literature keyword is not a measured binder; a patent association is not a structural complex. Preserve the weakest evidence category justified by the source.

## Target aliases

Alias changes in `config/target_aliases.json` must be conservative. Do not merge homologs, related family members, antigen domains, pathogen strains, or terms solely because they are commonly discussed together.

PLAbDab semicolon-delimited literature co-mentions must remain independent target interactions. A regression test in `tests/test_build_smoke.py` protects this rule.

`data/v4/target-review.json` is a curation aid only. Heavy antibody-set overlap can suggest an alias candidate, but it is not sufficient evidence for an automatic merge.

## Deduplication

PAIRS automatically merges exact normalized VH/VL pairs only. Do not introduce fuzzy sequence merging without a separately reviewed identity/provenance model.
