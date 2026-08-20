import csv
import json
from pathlib import Path

from pipeline.build import (
    _global_edit_metrics,
    _sequence_signature,
    compile_data,
    load_sources,
    merge_antibody,
    sequence_quality,
    write_indexes,
)
from pipeline.model import AntibodyObservation


def test_cluster_identity_boundaries_are_global_and_coverage_gated():
    reference = "ACDEFGHIKLMNPQRSTVWY" * 5
    at_90 = "Y" * 10 + reference[10:]
    below_90 = "Y" * 11 + reference[11:]
    assert _global_edit_metrics(reference, at_90) == (90.0, 100.0)
    assert _global_edit_metrics(reference, below_90) is None
    assert _global_edit_metrics(reference, reference[:89]) is None


def test_similarity_signature_retrieves_close_protein_candidates():
    reference = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWV"
    variant = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYVMSWVRQAPGKGLEWV"
    assert len(_sequence_signature(reference)) == 32
    assert set(_sequence_signature(reference)) & set(_sequence_signature(variant))


def test_sequence_quality_does_not_infer_completeness_from_length():
    quality = sequence_quality({"heavy": "EVQL" + "A" * 120, "light": "DIQM" + "A" * 105})
    assert quality["pairing"] == "paired"
    assert quality["completeness"] == "unknown_not_inferred"


def test_sequence_quality_requires_explicit_vhh_format():
    assert sequence_quality({"heavy": "AAA", "format": "VHH"})["explicit_vhh"] is True
    assert sequence_quality({"heavy": "AAA", "format": "antibody"})["explicit_vhh"] is False


def test_tiny_thera_build(tmp_path: Path):
    source = tmp_path / "thera.csv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Therapeutic",
                "HeavySequence",
                "LightSequence",
                "Target",
                "Format",
                "Est. Status",
                "Highest_Clin_Trial (Feb '25)",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "Therapeutic": "Examplemab",
                "HeavySequence": "AAA",
                "LightSequence": "BBB",
                "Target": "ERBB2/CD340/HER2",
                "Format": "Whole mAb",
                "Est. Status": "Active",
                "Highest_Clin_Trial (Feb '25)": "Phase-III",
            }
        )
    config = load_sources()
    stats = compile_data({"therasabdab": source}, config, tmp_path / "out")
    assert stats["antibodies"] == 1
    targets = json.loads((tmp_path / "out" / "targets.json").read_text())
    assert any(target["name"] == "ERBB2" for target in targets)
    assert (tmp_path / "out" / "sequence-search").exists()


def test_invalid_source_sequence_is_quarantined_without_target_claim(tmp_path: Path):
    source = tmp_path / "thera-invalid.csv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["Therapeutic", "HeavySequence", "LightSequence", "Target"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "Therapeutic": "Malformedmab",
                "HeavySequence": "AAA-CCC",
                "LightSequence": "BBB",
                "Target": "ERBB2",
            }
        )
    stats = compile_data({"therasabdab": source}, load_sources(), tmp_path / "out")
    assert stats["quarantined_records"] == 1
    assert stats["targets"] == 0
    antibody_file = next((tmp_path / "out" / "antibodies").glob("*.json"))
    antibody = next(iter(json.loads(antibody_file.read_text()).values()))
    assert antibody["metadata"]["sequence_quarantine"] is True
    assert antibody["metadata"]["sequence_quarantine_fields"]["heavy"]["raw"] == "AAA-CCC"


def test_source_nucleotide_stays_attached_to_exact_source_record():
    destination = {}
    merge_antibody(
        destination,
        AntibodyObservation(source="aa", record_id="record-a", name="same", heavy="AAA"),
    )
    merge_antibody(
        destination,
        AntibodyObservation(
            source="bb",
            record_id="record-b",
            name="same",
            heavy="AAA",
            vh_nt_source="ATGC",
            nucleotide_provenance={
                "VH": {
                    "source": "bb",
                    "source_record_id": "record-b",
                    "scope": "variable_domain",
                }
            },
            source_nucleotide_records=[
                {
                    "source": "bb",
                    "source_record_id": "record-b",
                    "chain": "VH",
                    "sequence": "ATGC",
                    "scope": "variable_domain",
                    "source_field": "vh_nt_source",
                }
            ],
        ),
    )
    records = {(
        record["source"],
        record["record_id"],
    ): record for record in destination["source_records"]}
    assert records[("aa", "record-a")]["nucleotide_records"] == []
    assert records[("bb", "record-b")]["nucleotide_records"][0]["sequence"] == "ATGC"


def test_plabdab_comentions_fan_out_without_alias_merging(tmp_path: Path):
    source = tmp_path / "plabdab.csv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["ID", "heavy_sequence", "light_sequence", "targets_mentioned"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "ID": "paper-antibody",
                "heavy_sequence": "EVQLVESGG",
                "light_sequence": "DIQMTQSPSS",
                "targets_mentioned": "ACP; ACP1; NDUFAB1; SDAP; FASN2A",
            }
        )
    config = load_sources()
    stats = compile_data({"plabdab": source}, config, tmp_path / "out")
    # Literature co-mentions are provenance/context only. They are retained
    # on the antibody record but cannot create primary target pages or direct
    # target assignments.
    assert stats["targets"] == 0
    assert stats["interactions"] == 5
    assert stats["literature_mentions"] == 5
    targets = json.loads((tmp_path / "out" / "targets.json").read_text())
    assert targets == []
    antibody_files = list((tmp_path / "out" / "antibodies").glob("*.json"))
    antibodies = json.loads(antibody_files[0].read_text())
    antibody = next(iter(antibodies.values()))
    assert antibody["direct_targets"] == []
    assert antibody["source_records"][0]["record_url"] == ""
    assert antibody["source_records"][0]["link_scope"] == "source_homepage"
    assert {item["target_name"] for item in antibody["literature_mentions"]} == {
        "ACP",
        "ACP1",
        "NDUFAB1",
        "SDAP",
        "FASN2A",
    }


def test_negative_only_observation_is_not_a_default_target_result(tmp_path: Path):
    stats = write_indexes(
        {
            "ab_negative": {
                "id": "ab_negative",
                "name": "Nonbinder",
                "sources": ["covabdab"],
                "aliases": [],
                "structures": [],
            }
        },
        [
            {
                "antibody_id": "ab_negative",
                "target_raw": "ERBB2",
                "relationship": "does_not_bind",
                "evidence": "CURATED",
                "source": "covabdab",
                "source_record_id": "Nonbinder",
                "reference": "",
                "epitope": "",
                "assay": "",
                "note": "",
            }
        ],
        tmp_path / "out",
        {"covabdab": {"records": 1, "interactions": 0}},
    )
    target = json.loads((tmp_path / "out" / "targets.json").read_text())[0]
    assert stats["targets"] == 1
    assert target["result_count"] == 0
    assert target["negative_count"] == 1
    target_dir = tmp_path / "out" / "targets" / target["dir"]
    assert not list(target_dir.glob("page-*.json"))
    negative_page = target_dir / "negative-page-001.json"
    assert json.loads(negative_page.read_text())[0]["relationships"] == ["does_not_bind"]


def test_functional_observation_is_separate_from_direct_target_results(tmp_path: Path):
    stats = write_indexes(
        {
            "ab_functional": {
                "id": "ab_functional",
                "name": "Neutralizer",
                "sources": ["covabdab"],
                "aliases": [],
                "structures": [],
            }
        },
        [
            {
                "antibody_id": "ab_functional",
                "target_raw": "SARS-CoV-2 Spike",
                "relationship": "neutralizes",
                "evidence": "CURATED",
                "source": "covabdab",
                "source_record_id": "Neutralizer",
                "reference": "",
                "epitope": "",
                "assay": "",
                "note": "",
            }
        ],
        tmp_path / "out",
        {"covabdab": {"records": 1, "interactions": 0}},
    )
    target = json.loads((tmp_path / "out" / "targets.json").read_text())[0]
    assert stats["targets"] == 1
    assert target["result_count"] == 0
    assert target["functional_count"] == 1
    assert target["page_count"] == 0
    assert target["functional_page_count"] == 1
    target_dir = tmp_path / "out" / "targets" / target["dir"]
    assert json.loads((target_dir / "functional-page-001.json").read_text())[0][
        "relationships"
    ] == ["neutralizes"]
