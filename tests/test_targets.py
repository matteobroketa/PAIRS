from pathlib import Path
import json

import pytest

from pipeline.targets import TargetResolver

ROOT = Path(__file__).resolve().parents[1]


def test_common_aliases_resolve_together():
    resolver = TargetResolver(ROOT / "config" / "target_aliases.json")
    assert resolver.resolve("HER2")[0] == resolver.resolve("ERBB2")[0]
    assert resolver.resolve("PD-L1")[0] == resolver.resolve("CD274")[0]
    assert resolver.resolve("SARS-CoV2 RBD")[0] == resolver.resolve("SARS-CoV-2 RBD")[0]
    assert resolver.resolve("Tumor necrosis factor")[0] == resolver.resolve("TNF-alpha")[0]


def test_muscle_achr_complex_and_subunits_keep_distinct_identities():
    resolver = TargetResolver(ROOT / "config" / "target_aliases.json")
    complex_id = resolver.resolve("AChR")[0]
    assert complex_id == resolver.resolve("acetylcholine receptor")[0]
    assert complex_id == resolver.resolve("muscle nAChR")[0]
    assert complex_id != resolver.resolve("CHRNA1")[0]
    assert resolver.resolve("acetylcholine receptor subunit alpha 1")[0] == resolver.resolve(
        "CHRNA1"
    )[0]
    assert resolver.entity("CHRNA1")["parent"] == "MUSCLE NICOTINIC ACETYLCHOLINE RECEPTOR"


def test_thera_synonym_group():
    resolver = TargetResolver(ROOT / "config" / "target_aliases.json")
    target_id, name, aliases = resolver.synonym_group("ERBB2/CD340/HER2")
    assert name == "ERBB2"
    assert "HER2" in aliases
    assert target_id == resolver.resolve("ERBB2")[0]


def test_literature_comentions_never_become_aliases():
    resolver = TargetResolver(ROOT / "config" / "target_aliases.json")
    groups = resolver.mention_group("ACP; ACP1; NDUFAB1; SDAP; FASN2A")
    assert len(groups) == 5
    ids = [group[0] for group in groups]
    assert len(set(ids)) == 5
    for _, _, aliases in groups:
        assert (
            len(
                {
                    alias
                    for alias in aliases
                    if alias in {"ACP", "ACP1", "NDUFAB1", "SDAP", "FASN2A"}
                }
            )
            == 1
        )


def test_unresolved_target_ids_do_not_collide_when_slugs_match():
    resolver = TargetResolver(ROOT / "config" / "target_aliases.json")
    assert resolver.resolve("A/B")[0] != resolver.resolve("A B")[0]


def test_authoritative_mapping_does_not_change_local_target_identity():
    resolver = TargetResolver(ROOT / "config" / "target_aliases.json")
    target_id, canonical, _ = resolver.resolve("HER2")
    entity = resolver.entity(canonical)
    assert target_id == "target:erbb2"
    assert entity["identifiers"]["uniprot"] == "P04626"
    assert entity["mapping_provenance"]["scope"] == "exact canonical target only"


def test_stable_identifiers_resolve_before_curated_free_text_aliases():
    resolver = TargetResolver(ROOT / "config" / "target_aliases.json")
    assert resolver.resolve("unrelated antigen", "UniProt:P00533")[0] == resolver.resolve("EGFR")[0]
    assert resolver.resolve(
        "unrelated antigen", "https://www.uniprot.org/uniprotkb/P00533/entry"
    )[0] == resolver.resolve("EGFR")[0]
    assert resolver.resolve("unrelated antigen", "HGNC:3430")[0] == resolver.resolve("ERBB2")[0]
    assert resolver.resolve("unrelated antigen", "NCBI Gene:1134")[0] == resolver.resolve(
        "CHRNA1"
    )[0]


def test_name_resolution_keeps_species_mismatch_unresolved_but_id_wins():
    resolver = TargetResolver(ROOT / "config" / "target_aliases.json")
    unresolved = resolver.resolve("EGFR", organism="Mus musculus")
    assert unresolved[0] != resolver.resolve("EGFR")[0]
    assert "Mus musculus" in unresolved[1]
    assert resolver.resolve("EGFR", "UniProt:P00533", "Mus musculus")[0] == resolver.resolve(
        "EGFR"
    )[0]


def test_unknown_stable_identifier_is_retained_as_non_colliding_identity():
    resolver = TargetResolver(ROOT / "config" / "target_aliases.json")
    first = resolver.resolve("first display label", "UniProt:Q9TEST1")[0]
    second = resolver.resolve("different display label", "UniProt:Q9TEST1")[0]
    assert first == second
    assert first.startswith("target:identifier-uniprot-q9test1")


def test_target_alias_collision_fails_closed(tmp_path):
    aliases = tmp_path / "aliases.json"
    aliases.write_text(json.dumps({"Target A": ["shared"], "Target B": ["SHARED"]}))
    with pytest.raises(ValueError, match="target alias collision"):
        TargetResolver(aliases)
