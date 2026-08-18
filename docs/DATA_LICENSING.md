# Data licensing and redistribution policy

The repository's original **software** is MIT licensed. Imported datasets are separate works and retain their upstream licenses/terms.

The build pipeline is source-aware. `config/sources.json` records the upstream homepage, citation, licensing note, and build-time discovery rule. Restricted or permission-based resources are disabled in the default public build rather than copied into the repository.

Before publishing a large public mirror, maintainers should independently confirm the current redistribution terms for every enabled source. Upstream terms can change. When in doubt, ship an adapter that users run against their own legitimately obtained copy instead of redistributing the dataset.

The tracked `data/v4/` directory is a compact integration/starter snapshot and should not be treated as a canonical archival copy of any upstream database. A scheduled or manually requested GitHub Actions workflow builds and validates the uncapped public snapshot. Normal code deployments reuse its latest successful production-data artifact and never rebuild the corpus.
