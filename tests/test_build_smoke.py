import csv
import json
from pathlib import Path

import pipeline.build as build_module
import pipeline.sources as source_module
from pipeline.build import (
    _global_edit_metrics,
    _parse_accession_sequence,
    _observation_accession_keys,
    _sequence_signature,
    compile_data,
    download_nanosaurus_json,
    load_sources,
    merge_antibody,
    recover_accession_records,
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


def test_same_stable_identifier_allows_source_display_label_variants(tmp_path: Path):
    antibodies = {
        "ab_a": {"id": "ab_a", "name": "A", "heavy": "AAA", "light": "BBB"},
        "ab_b": {"id": "ab_b", "name": "B", "heavy": "CCC", "light": "DDD"},
    }
    interactions = [
        {
            "antibody_id": "ab_a",
            "target_raw": "NLRP1",
            "target_external_id": "UniProt:Q9C000",
            "relationship": "binds",
            "evidence": "CURATED",
            "source": "source_a",
            "source_record_id": "a1",
        },
        {
            "antibody_id": "ab_b",
            "target_raw": "NACHT, LRR and PYD domains-containing protein 1",
            "target_external_id": "UniProt:Q9C000",
            "relationship": "binds",
            "evidence": "CURATED",
            "source": "source_b",
            "source_record_id": "b1",
        },
    ]
    stats = write_indexes(
        antibodies,
        interactions,
        tmp_path / "out",
        {"source_a": {"records": 1}, "source_b": {"records": 1}},
    )
    assert stats["targets"] == 1
    targets = json.loads((tmp_path / "out" / "targets.json").read_text())
    assert len(targets[0]["aliases"]) >= 2


def test_nanosaurus_downloader_paginates_public_json(monkeypatch, tmp_path: Path):
    pages = {
        1: {"response_type": "success", "page": 1, "max_pages": 2, "results": [{"id": 1}]},
        2: {"response_type": "success", "page": 2, "max_pages": 2, "results": [{"id": 2}]},
    }

    class Response:
        def __init__(self, payload):
            self.payload = payload

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def read(self):
            return json.dumps(self.payload).encode("utf-8")

    def fake_urlopen(request, timeout):
        page = json.loads(request.data.decode("utf-8"))["page"]
        return Response(pages[page])

    monkeypatch.setattr(build_module.urllib.request, "urlopen", fake_urlopen)
    path = tmp_path / "nanosaurus.json"
    info = download_nanosaurus_json("https://example.test/api/query", path, attempts=1)
    assert info["pages"] == 2
    assert info["records"] == 2
    assert json.loads(path.read_text())["results"] == [{"id": 1}, {"id": 2}]


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


def test_targeted_accession_recovery_caches_exact_record_and_provenance(tmp_path):
    observation = AntibodyObservation(
        source="curated",
        record_id="record-1",
        name="Example",
        sequence_accessions=[
            {
                "accession": "X77017.1",
                "database": "genbank",
                "chain": "VH",
                "scope": "variable_domain",
                "source_field": "VH GenBank accession",
            }
        ],
    )
    payload = "LOCUS       X77017       12 bp\nORIGIN\n        1 atgc ccnn ataa\n//\n"
    assert recover_accession_records(observation, tmp_path, lambda url: payload) == 1
    assert observation.vh_nt_source == "ATGCCCNNATAA"
    assert observation.sequence_architecture == "nucleotide_only"
    assert observation.metadata["sequence_architecture"] == "nucleotide_only"
    assert observation.nucleotide_provenance["VH"]["source_record_id"] == "X77017.1"
    assert observation.source_nucleotide_records[0]["source"] == "genbank_accession"
    assert list((tmp_path / "accessions" / "genbank").glob("X77017.1.json"))

    def unexpected_fetch(url):
        raise AssertionError("cached accession was fetched again")

    assert recover_accession_records(observation, tmp_path, unexpected_fetch) == 1


def test_accession_parser_accepts_genbank_ena_fasta_and_embl_formats():
    assert _parse_accession_sequence(
        "genbank", "LOCUS x\nORIGIN\n        1 atg ccnn\n//\n"
    ) == "ATGCCNN"
    assert _parse_accession_sequence("ena", ">X1\natg ccnn\n") == "ATGCCNN"
    assert _parse_accession_sequence(
        "embl", "ID   X1; SV 1;\nSQ   Sequence 7 BP;\n     atg ccnn\n//\n"
    ) == "ATGCCNN"


def test_cross_source_reconciliation_does_not_merge_heavy_only_chains(monkeypatch, tmp_path):
    def adapter_a(path, source_url=""):
        yield AntibodyObservation(
            source="source_a",
            record_id="a1",
            name="A",
            heavy="EVQLV",
            sequence_accessions=[{"database": "genbank", "accession": "X1", "chain": "VH"}],
        ), []

    def adapter_b(path, source_url=""):
        yield AntibodyObservation(
            source="source_b",
            record_id="b1",
            name="B",
            heavy="EVQLV",
            sequence_accessions=[{"database": "embl", "accession": "X1", "chain": "VH"}],
        ), []

    monkeypatch.setitem(source_module.ADAPTERS, "test_a", adapter_a)
    monkeypatch.setitem(source_module.ADAPTERS, "test_b", adapter_b)
    config = {
        "source_a": {"adapter": "test_a"},
        "source_b": {"adapter": "test_b"},
    }
    stats = compile_data(
        {"source_a": tmp_path / "a", "source_b": tmp_path / "b"},
        config,
        tmp_path / "out",
    )
    assert stats["antibodies"] == 2
    assert stats["reconciliation_counts"] == {}


def test_insdc_accession_namespace_is_database_independent():
    genbank = AntibodyObservation(
        source="a",
        record_id="1",
        name="A",
        sequence_accessions=[{"database": "genbank", "accession": "X1", "chain": "VH"}],
    )
    embl = AntibodyObservation(
        source="b",
        record_id="2",
        name="B",
        sequence_accessions=[{"database": "embl", "accession": "X1", "chain": "VH"}],
    )
    assert _observation_accession_keys(genbank) == _observation_accession_keys(embl)


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


def test_plabdab_sequence_survives_empty_target_metadata(tmp_path: Path):
    source = tmp_path / "plabdab-empty-target.csv"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["ID", "heavy_sequence", "light_sequence", "targets_mentioned"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "ID": "sequence-without-target",
                "heavy_sequence": "EVQLVESGG",
                "light_sequence": "DIQMTQSPSS",
                "targets_mentioned": "",
            }
        )

    stats = compile_data({"plabdab": source}, load_sources(), tmp_path / "out")
    assert stats["antibodies"] == 1
    assert stats["interactions"] == 0
    assert stats["targets"] == 0
    shard = next((tmp_path / "out" / "antibodies").glob("*.json"))
    antibody = next(iter(json.loads(shard.read_text()).values()))
    assert antibody["heavy"] == "EVQLVESGG"
    assert antibody["light"] == "DIQMTQSPSS"


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
