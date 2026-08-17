import csv
import json
from pathlib import Path

from pipeline.build import compile_data, load_sources
from pipeline.sources import therasabdab


THERA_FIELDS = [
    "Therapeutic",
    "HeavySequence",
    "LightSequence",
    "HeavySequence(ifbispec)",
    "LightSequence(ifbispec)",
    "Target",
    "Format",
]


def _write_thera(path: Path, row: dict) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=THERA_FIELDS)
        writer.writeheader()
        writer.writerow(row)


def test_thera_multispecific_preserves_arms_without_inferred_target_claims(tmp_path: Path):
    source = tmp_path / "thera.csv"
    _write_thera(
        source,
        {
            "Therapeutic": "Acapatamab",
            "HeavySequence": "AAA",
            "LightSequence": "BBB",
            "HeavySequence(ifbispec)": "CCC",
            "LightSequence(ifbispec)": "DDD",
            "Target": "CD3E;BCMA",
            "Format": "Bispecific mAb",
        },
    )

    observations = list(therasabdab(source))
    assert len(observations) == 2
    assert all(not interactions for _, interactions in observations)
    assert {antibody.arm["label"] for antibody, _ in observations} == {"arm_1", "arm_2"}
    assert all(antibody.metadata["multispecific"] is True for antibody, _ in observations)
    assert all(
        antibody.metadata["target_assignment_status"] == "unavailable_no_arm_mapping"
        for antibody, _ in observations
    )
    assert observations[0][0].construct["target_raw"] == ["CD3E", "BCMA"]
    assert all("target_raw" not in antibody.arm for antibody, _ in observations)

    stats = compile_data({"therasabdab": source}, load_sources(), tmp_path / "out")
    assert stats["antibodies"] == 2
    assert stats["source_counts"]["therasabdab"]["records"] == 1
    assert stats["interactions"] == 0
    assert stats["targets"] == 0
    assert json.loads((tmp_path / "out" / "targets.json").read_text()) == []
    shard_records = []
    for shard in (tmp_path / "out" / "antibodies").glob("*.json"):
        shard_records.extend(json.loads(shard.read_text()).values())
    assert all(record["constructs"][0]["target_assignment_status"] == "unavailable_no_arm_mapping" for record in shard_records)
    assert all(record.get("therapeutic_status", "") == "" and record.get("structures", []) == [] for record in shard_records)


def test_thera_ambiguous_target_delimiter_is_quarantined_without_format_marker(tmp_path: Path):
    source = tmp_path / "thera.csv"
    _write_thera(
        source,
        {
            "Therapeutic": "Example-ambiguous",
            "HeavySequence": "AAA",
            "LightSequence": "BBB",
            "Target": "CD3E;BCMA",
            "Format": "Whole mAb",
        },
    )

    antibody, interactions = list(therasabdab(source))[0]
    assert interactions == []
    assert antibody.metadata["target_assignment_status"] == "unavailable_no_arm_mapping"


def test_thera_structure_identity_tiers_keep_exact_separate_from_homologs(tmp_path: Path):
    source = tmp_path / "thera.csv"
    fields = THERA_FIELDS + ["100% SI Structure", "99% SI Structure", "95-98% SI Structure"]
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow(
            {
                "Therapeutic": "Tieredmab",
                "HeavySequence": "AAA",
                "LightSequence": "BBB",
                "Target": "ERBB2/CD340/HER2",
                "Format": "Whole mAb",
                "100% SI Structure": "8QYB:HL",
                "99% SI Structure": "8QYA:HL/8QY9:HL",
                "95-98% SI Structure": "8QY8:HL",
            }
        )

    antibody, interactions = list(therasabdab(source))[0]
    assert interactions[0].target_raw == "ERBB2/CD340/HER2"
    assert antibody.structures == ["8QYB"]
    assert antibody.structure_tiers == {
        "100%": ["8QYB"],
        "99%": ["8QY9", "8QYA"],
        "95-98%": ["8QY8"],
    }

    compile_data({"therasabdab": source}, load_sources(), tmp_path / "out")
    shard = next((tmp_path / "out" / "antibodies").glob("*.json"))
    record = next(iter(json.loads(shard.read_text()).values()))
    assert record["structures"] == ["8QYB"]
    assert record["structure_tiers"]["99%"] == ["8QY9", "8QYA"]


def test_thera_provenance_uses_authoritative_record_summary_url(tmp_path: Path):
    source = tmp_path / "thera.csv"
    _write_thera(
        source,
        {
            "Therapeutic": "Examplemab",
            "HeavySequence": "AAA",
            "LightSequence": "BBB",
            "Target": "ERBB2",
            "Format": "Whole mAb",
        },
    )
    antibody, interactions = list(
        therasabdab(source, "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/therasabdab/search/")
    )[0]
    expected = (
        "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/therasabdab/therasummary/?INN=Examplemab"
    )
    assert antibody.record_url == expected
    assert antibody.link_scope == "record"
    assert interactions[0].record_url == ""

    compile_data({"therasabdab": source}, load_sources(), tmp_path / "out")
    shard = next((tmp_path / "out" / "antibodies").glob("*.json"))
    record = next(iter(json.loads(shard.read_text()).values()))
    assert record["source_records"][0]["record_url"] == expected
    assert record["source_records"][0]["link_scope"] == "record"
    target = json.loads((tmp_path / "out" / "targets.json").read_text())[0]
    page = tmp_path / "out" / "targets" / target["dir"] / "page-001.json"
    row = json.loads(page.read_text())[0]
    assert row["interactions"][0]["record_url"] == expected
    assert row["interactions"][0]["link_scope"] == "record"
