import json
from pathlib import Path

from pipeline.build import compile_data, load_sources
from pipeline.sources import _sabdab_targets, sabdab2


def _snapshot(path: Path) -> None:
    path.write_text(
        json.dumps(
            {
                "records": [
                    {
                        "id": "sabdab2_H0001L0001",
                        "heavy_like_segment": {
                            "sequence": "EVQLVESGG",
                            "cdr1": "GFTF",
                            "cdr2": "ISW",
                            "cdr3": "ARDR",
                        },
                        "light_like_segment": {
                            "sequence": "DIQMTQSPSS",
                            "cdr1": "SSN",
                            "cdr2": "YAS",
                            "cdr3": "QQY",
                        },
                        "instances": [{"id": "pdb_00005hbt-A-B", "type": "FAB"}],
                    },
                    {
                        "id": "sabdab2_EGFR",
                        "heavy_like_segment": {"sequence": "EVQLVAAEGFR"},
                        "light_like_segment": {"sequence": "DIQMTAAEGFR"},
                        "instances": [{"id": "pdb_00006egf-A-B", "type": "FAB"}],
                    },
                    {
                        "id": "sabdab2_TNF",
                        "heavy_like_segment": {"sequence": "EVQLVAATNF"},
                        "light_like_segment": {"sequence": "DIQMT AATNF"},
                        "instances": [{"id": "pdb_00007tnf-A-B", "type": "FAB"}],
                    },
                ],
                "summary": [
                    {
                        "INSTANCE": "pdb_00005hbt-A-B",
                        "PDB": "pdb_00005hbt",
                        "SABDAB_ID": "sabdab2_H0001L0001",
                        "antigen_chain": "B",
                        "antigen_type": "PROTEIN",
                        "antigen_name": "Acetylcholine receptor subunit alpha 1",
                        "compound": "Fab35 structural study",
                        "organism": "rattus norvegicus",
                        "type": "FAB",
                        "antigen_uniprot": "P02708",
                    },
                    {
                        "INSTANCE": "pdb_00006egf-A-B",
                        "PDB": "pdb_00006egf",
                        "SABDAB_ID": "sabdab2_EGFR",
                        "antigen_chain": "B",
                        "antigen_type": "PROTEIN",
                        "antigen_name": "epidermal growth factor receptor",
                        "antigen_uniprot": "P00533",
                        "compound": "FabEGFR structural study",
                        "organism": "Homo sapiens",
                        "type": "FAB",
                    },
                    {
                        "INSTANCE": "pdb_00007tnf-A-B",
                        "PDB": "pdb_00007tnf",
                        "SABDAB_ID": "sabdab2_TNF",
                        "antigen_chain": "B",
                        "antigen_type": "PROTEIN",
                        "antigen_name": "tumor necrosis factor",
                        "antigen_uniprot": "P01375",
                        "compound": "FabTNF structural study",
                        "organism": "Homo sapiens",
                        "type": "FAB",
                    },
                ],
            }
        ),
        encoding="utf-8",
    )


def test_sabdab2_maps_sequences_and_explicit_antigen_pairings(tmp_path: Path):
    source = tmp_path / "sabdab2.json"
    _snapshot(source)

    antibody, interactions = next(sabdab2(source))
    assert antibody.heavy == "EVQLVESGG"
    assert antibody.light == "DIQMTQSPSS"
    assert antibody.name == "Fab35"
    assert "Fab35" in antibody.aliases
    assert antibody.structures == ["5HBT"]
    assert antibody.record_url.endswith("/pdb/5hbt/")
    assert len(interactions) == 1
    assert interactions[0].relationship == "binds"
    assert interactions[0].target_raw == "Acetylcholine receptor subunit alpha 1"
    assert interactions[0].evidence == "STRUCTURE"

    stats = compile_data({"sabdab2": source}, load_sources(), tmp_path / "out")
    assert stats["antibodies"] == 3
    assert stats["interactions"] == 3
    assert stats["targets"] == 3
    targets = json.loads((tmp_path / "out" / "targets.json").read_text())
    by_name = {target["name"]: target for target in targets}
    assert {"CHRNA1", "EGFR", "TNF"} == set(by_name)
    assert by_name["CHRNA1"]["external_ids"] == ["UniProt:P02708"]
    assert "P02708" in by_name["CHRNA1"]["aliases"]
    assert "UNIPROT:P02708" in by_name["CHRNA1"]["aliases"]
    assert "HGNC:1955" in by_name["CHRNA1"]["aliases"]
    assert "NCBI_GENE:1134" in by_name["CHRNA1"]["aliases"]
    assert by_name["EGFR"]["external_ids"] == ["UniProt:P00533"]
    assert by_name["TNF"]["external_ids"] == ["UniProt:P01375"]


def test_sabdab2_does_not_broadcast_single_identifier_across_antigens():
    targets = _sabdab_targets(
        {
            "antigen_name": "Protein A|Protein B",
            "antigen_uniprot": "P12345",
            "antigen_species": "Homo sapiens",
        }
    )
    assert [target["external_id"] for target in targets] == ["", ""]
    assert [target["organism"] for target in targets] == ["Homo sapiens", "Homo sapiens"]


def test_sabdab2_keeps_positional_identifiers_aligned():
    targets = _sabdab_targets(
        {
            "antigen_name": "Protein A|Protein B",
            "antigen_uniprot": "P12345|P67890",
        }
    )
    assert [target["external_id"] for target in targets] == [
        "UniProt:P12345",
        "UniProt:P67890",
    ]
