# Similarity benchmarks

PAIRS reports candidate-stage recall separately from alignment accuracy. Run:

```bash
python -m pipeline.benchmark data/v4 --samples 100 --output similarity-benchmark.json
```

The benchmark selects a deterministic sequence sample, creates controlled substitution relatives at
100%, 99%, 95%, 90%, and 80% identity, and records whether the indexed 5-mer stage retrieves the
known source sequence. Exact queries include the production SHA-256 fallback; mutated queries
measure candidate recall before alignment. This is not a biological lineage benchmark and does not
establish shared specificity. Recall below 100% is an explicit false-negative estimate, not a passing
score hidden by alignment.

Controlled variants are a reproducible baseline. A separately curated set of natural antibody
relatives remains necessary for external biological validation.
