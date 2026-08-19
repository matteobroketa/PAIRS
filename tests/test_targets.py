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


def test_target_alias_collision_fails_closed(tmp_path):
    aliases = tmp_path / "aliases.json"
    aliases.write_text(json.dumps({"Target A": ["shared"], "Target B": ["SHARED"]}))
    with pytest.raises(ValueError, match="target alias collision"):
        TargetResolver(aliases)
