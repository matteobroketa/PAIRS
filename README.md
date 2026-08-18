# PAIRS

**Pan-Antibody Integrated Retrieval System**

Search public antibody records by target, antibody name or sequence.

[Open PAIRS](https://matteobroketa.github.io/PAIRS/)

PAIRS brings antibody sequences, targets, structures and source evidence into one searchable interface. Results remain linked to their original databases so important records can be checked at the source.

## What you can do

- Find antibodies reported against a target
- Search antibody and therapeutic names
- Search VH/VHH, VL and CDR sequences
- Retrieve full-chain VH/VL similarity candidates with explicit alignments
- Screen multiple sequences against the public corpus
- Browse targets with therapeutic, structural or paired-sequence evidence
- Inspect antibody sequences, linked measurements, source annotations, and lazy embedded PDB views
- Compare 2–10 selected antibodies with explicit global amino-acid alignments
- Build source-record, exact-pair, or 99/95/90 sequence-cluster-deduplicated datasets
- Filter by evidence source, reported origin and sequence availability
- Save selections and recent text views in a local-only workspace
- Screen multi-FASTA or paired `name,VH,VL` CSV/TSV panels locally
- Export results as FASTA, metadata/provenance, and BibTeX/RIS/CSV citations
- Open stable target/antibody entity views with explicit hierarchy and provenance timelines
- Query the static snapshot from Python environments with the local `pairs` CLI
- Generate clearly separated, parameterized synthetic coding-DNA bundles

## Data sources

PAIRS currently indexes public records from:

- PLAbDab
- Thera-SAbDab
- CoV-AbDab
- Pox-AbDab
- IEDB BCR

Each result retains its source and evidence type. See [docs/SOURCES.md](docs/SOURCES.md) for source details, citations and data-use terms.

## Data note

PAIRS is a search and indexing tool, not an experimental validation database. Records may contain incomplete sequences or source-specific annotations; verify important findings against the linked original source.

Sequence similarity results are retrieval candidates and do not by themselves establish antibody identity or shared antigen specificity.
Sequence clusters are indexed redundancy aids, not clonal lineages.

## Run locally

PAIRS is a static site.

```bash
python -m http.server 8000
```

Then open:

```text
http://localhost:8000
```

## Rebuild the data

Install the Python package and rebuild the static search data:

```bash
pip install -e .
python -m pipeline.build
python -m pipeline.validate
```

The included GitHub Actions workflows can rebuild the data and deploy the site to GitHub Pages.

## Programmatic data access

PAIRS exposes its search data as static JSON files that can also be used from scripts and other tools.

See [docs/DATA_ACCESS.md](docs/DATA_ACCESS.md).

## Contributing

Bug reports, source corrections and contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md).

## Citation

Citation metadata are provided in [CITATION.cff](CITATION.cff).

## License

See [LICENSE](LICENSE) for the PAIRS software license. Individual source datasets retain their own terms and attribution requirements.
