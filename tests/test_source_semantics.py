import csv
import json

from pipeline.build import compile_data, load_sources
from pipeline.sources import _cov_target, _pox_entity, iedb


def test_cov_parent_target_is_labelled_as_derived_hierarchy():
    targets = dict(_cov_target("SARS-CoV-2", "RBD"))
    assert targets["SARS-CoV-2 RBD"] == "source_epitope"
    assert targets["SARS-CoV-2 Spike"] == "derived_hierarchy"


def test_pox_target_retains_strain_location_and_protein():
    target, raw = _pox_entity("VACV_WR_Copenhagen_A27")
    assert target == "Vaccinia virus WR Copenhagen A27"
    assert raw == "VACV_WR_Copenhagen_A27"


def test_iedb_keeps_chain_annotation_provenance_and_linkage(tmp_path):
    path = tmp_path / "iedb.csv"
    row = {
        "receptor_group_id": "2191",
        "receptor__group_iri": "https://www.iedb.org/receptor/2191",
        "receptor__iedb_receptor_id": "IEDB_RECEPTOR:1",
        "receptor__reference_name": "Example receptor",
        "receptor__type": "paired antibody",
        "reference__iedb_iri": "https://www.iedb.org/reference/12",
        "epitope__iedb_iri": "https://www.iedb.org/epitope/34",
        "epitope__name": "peptide ABC",
        "epitope__source_molecule": "Example protein",
        "assay__type": "binding assay",
        "assay__iedb_ids": "IEDB:55;IEDB:56",
        "chain_1__type": "heavy",
        "chain_1__nucleotide_sequence": "ATGCNN",
        "chain_1__protein_sequence": "FULLHEAVY",
        "chain_1__v_domain_calculated": "EVQLVESGG",
        "chain_1__cdr1_curated": "OLD",
        "chain_1__cdr1_calculated": "CALC",
        "chain_1__cdr2_calculated": "WINT",
        "chain_1__cdr3_calculated": "CAR",
        "chain_1__calculated_v_gene": "IGHV1",
        "chain_1__calculated_d_gene": "IGHD1",
        "chain_1__calculated_j_gene": "IGHJ1",
        "chain_2__type": "lambda light",
        "chain_2__v_domain_calculated": "DIQMTQ",
        "chain_2__cdr3_curated": "CQY",
        "chain_2__curated_v_gene": "IGLV1",
        "chain_2__curated_j_gene": "IGLJ1",
    }
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)

    antibody, interactions = next(iedb(path, "https://www.iedb.org/database_export_v3.php"))
    assert antibody.heavy == "EVQLVESGG"
    assert antibody.light == "DIQMTQ"
    assert antibody.cdrh1 == "CALC"
    assert antibody.heavy_d == "IGHD1"
    assert antibody.chain_annotations["VH"]["source_full_protein"] == "FULLHEAVY"
    assert antibody.chain_annotations["VH"]["regions"]["CDR1"]["curated"] == "OLD"
    assert antibody.nucleotide_provenance["VH"]["scope"] == "full_length_bcr"
    # An epitope IRI must not be mislabeled as the source-molecule identifier.
    assert interactions[0].target_external_id == ""
    assert interactions[0].assay_ids == ["IEDB:55", "IEDB:56"]
    assert interactions[0].receptor_group_id == "2191"

    stats = compile_data({"iedb": path}, load_sources(), tmp_path / "out")
    assert stats["interactions"] == 1
    target = json.loads((tmp_path / "out" / "targets.json").read_text())[0]
    page = json.loads((tmp_path / "out" / "targets" / target["dir"] / "page-001.json").read_text())
    built = page[0]["interactions"][0]
    assert built["target_external_id"] == ""
    assert built["assay_ids"] == ["IEDB:55", "IEDB:56"]


def test_iedb_sequence_only_record_does_not_invent_target(tmp_path):
    path = tmp_path / "iedb-empty-target.csv"
    row = {
        "receptor_group_id": "7",
        "receptor__type": "BCR",
        "epitope__name": "",
        "chain_1__type": "heavy",
        "chain_1__v_domain_calculated": "EVQLV",
        "chain_2__type": "",
    }
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
    antibody, interactions = next(iedb(path))
    assert antibody.heavy == "EVQLV"
    assert antibody.light == ""
    assert interactions == []


def test_iedb_measurements_join_only_by_explicit_group_and_assay(tmp_path):
    receptor_path = tmp_path / "iedb.csv"
    receptor = {
        "receptor_group_id": "91",
        "receptor__type": "heavylight",
        "epitope__name": "",
        "chain_1__type": "heavy",
        "chain_1__v_domain_calculated": "EVQLVESGG",
        "chain_2__type": "kappa_light",
        "chain_2__v_domain_calculated": "DIQMTQ",
    }
    with receptor_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(receptor))
        writer.writeheader()
        writer.writerow(receptor)

    search_path = tmp_path / "iedb-bcr-search.csv"
    with search_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["receptor_group_id", "iedb_assay_ids"])
        writer.writeheader()
        writer.writerow({"receptor_group_id": "91", "iedb_assay_ids": "{501,502}"})
        writer.writerow({"receptor_group_id": "92", "iedb_assay_ids": "{999}"})

    measurement_path = tmp_path / "iedb-bcell-export.csv"
    fields = [
        "assay_id",
        "assay_id__iedb_iri",
        "reference__iedb_iri",
        "epitope__name",
        "epitope__molecule_parent",
        "epitope__molecule_parent_iri",
        "assay__method",
        "assay__response_measured",
        "assay__units",
        "assay__qualitative_measure",
        "assay__measurement_inequality",
        "assay__quantitative_measurement",
        "assay__comments",
    ]
    with measurement_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow(
            {
                "assay_id": "501",
                "assay_id__iedb_iri": "https://www.iedb.org/assay/501",
                "epitope__name": "epitope A",
                "epitope__molecule_parent": "Example antigen",
                "epitope__molecule_parent_iri": "UNIPROT:P00001",
                "assay__method": "surface plasmon resonance",
                "assay__response_measured": "dissociation constant KD",
                "assay__units": "nM",
                "assay__qualitative_measure": "Positive",
                "assay__measurement_inequality": "<",
                "assay__quantitative_measurement": "2.5",
            }
        )
        writer.writerow(
            {
                "assay_id": "502",
                "assay_id__iedb_iri": "https://www.iedb.org/assay/502",
                "epitope__name": "epitope A",
                "epitope__molecule_parent": "Example antigen",
                "epitope__molecule_parent_iri": "UNIPROT:P00001",
                "assay__method": "ELISA",
                "assay__response_measured": "half maximal effective concentration EC50",
                "assay__units": "nM",
                "assay__qualitative_measure": "Positive",
                "assay__quantitative_measurement": "9",
            }
        )
        writer.writerow(
            {
                "assay_id": "999",
                "epitope__molecule_parent": "Wrong receptor antigen",
                "assay__response_measured": "dissociation constant KD",
                "assay__units": "nM",
                "assay__quantitative_measurement": "1",
            }
        )

    _, interactions = next(iedb(receptor_path))
    assert [item.assay_ids for item in interactions] == [["501"], ["502"]]
    assert [item.measurements[0]["metric"] for item in interactions] == ["KD", "EC50"]
    assert interactions[0].measurements[0]["qualifier"] == "<"
    assert {item.target_raw for item in interactions} == {"Example antigen"}
