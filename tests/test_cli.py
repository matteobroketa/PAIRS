import json

from pipeline.build import write_indexes
from pipeline.cli import main


def _snapshot(tmp_path):
    output = tmp_path / "v4"
    antibodies = {
        "ab_parent": {
            "id": "ab_parent",
            "name": "Spike antibody",
            "heavy": "EVQLVESGG",
            "light": "DIQMTQ",
            "sources": ["covabdab"],
            "aliases": [],
            "structures": [],
        },
        "ab_child": {
            "id": "ab_child",
            "name": "RBD antibody",
            "heavy": "QVQLVQSGA",
            "light": "EIVLTQSPA",
            "sources": ["covabdab"],
            "aliases": [],
            "structures": [],
        },
    }
    interactions = [
        {
            "antibody_id": "ab_parent",
            "target_raw": "SARS-CoV-2 Spike",
            "relationship": "binds",
            "evidence": "CURATED",
            "source": "covabdab",
            "source_record_id": "parent",
            "reference": "",
            "epitope": "",
            "assay": "",
            "note": "",
        },
        {
            "antibody_id": "ab_child",
            "target_raw": "SARS-CoV-2 RBD",
            "relationship": "binds",
            "evidence": "CURATED",
            "source": "covabdab",
            "source_record_id": "child",
            "reference": "",
            "epitope": "",
            "assay": "",
            "note": "",
        },
    ]
    stats = write_indexes(
        antibodies,
        interactions,
        output,
        {"covabdab": {"records": 2, "interactions": 0}},
    )
    (output / "manifest.json").write_text(
        json.dumps({"schema_version": 4, "snapshot": "test", "stats": stats}),
        encoding="utf-8",
    )
    return output


def test_cli_target_hierarchy_is_explicit(tmp_path, capsys):
    snapshot = _snapshot(tmp_path)
    assert main(["--data", str(snapshot), "target", "SARS-CoV-2 Spike"]) == 0
    exact = json.loads(capsys.readouterr().out)
    assert exact["count"] == 1
    assert exact["query"]["scope"] == "exact"

    assert (
        main(
            [
                "--data",
                str(snapshot),
                "target",
                "SARS-CoV-2 Spike",
                "--include-descendants",
            ]
        )
        == 0
    )
    expanded = json.loads(capsys.readouterr().out)
    assert expanded["count"] == 2
    assert expanded["query"]["scope"] == "include_descendants"


def test_cli_exact_sequence_and_fasta(tmp_path, capsys):
    snapshot = _snapshot(tmp_path)
    assert (
        main(
            [
                "--data",
                str(snapshot),
                "--format",
                "fasta",
                "sequence",
                "EVQLVESGG",
            ]
        )
        == 0
    )
    output = capsys.readouterr().out
    assert "pairs_id=ab_parent|chain=VH|sequence=amino_acid" in output


def test_cli_stale_id_fails_without_substitution(tmp_path, capsys):
    snapshot = _snapshot(tmp_path)
    assert main(["--data", str(snapshot), "antibody", "ab_missing"]) == 2
    assert "not present in this snapshot" in capsys.readouterr().err
