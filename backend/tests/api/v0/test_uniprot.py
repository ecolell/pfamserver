import json


def test_get_uniprot_pfams(db, client, uniprot_reg_full_mt2_human):
    headers = [("Accept", "application/json"), ("Content-Type", "application/json")]
    res = client.get("/api/v0/uniprots/mt2_human/pfams", headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data["query"] == "MT2_HUMAN"
    assert data["output"] == [
        {
            "num_full": "480",
            "description": "Metallothionein",
            "seq_start": 1,
            "pfamA_acc": "PF00131",
            "seq_end": 61,
        }
    ]

    res = client.get("/api/v0/uniprots/invalid_uniprot/pfams", headers=headers)
    assert res.status_code == 400
    data = json.loads(res.get_data(as_text=True))
    assert data["message"] == "Uniprot doesn" "t exist."


def test_get_reference_sequences_pdbs(
    db, client, uniprot_mt2_human, pdb_pfam_a_reg_pf01030
):
    headers = [("Accept", "application/json"), ("Content-Type", "application/json")]
    res = client.get("/api/v0/uniprots/MT2_HUMAN/1-61/pdbs", headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data["query"] == {"seq_end": 61, "seq_start": 1, "uniprot_id": "MT2_HUMAN"}
    assert len(data["output"]) == 2


def test_get_uniprots(db, client, uniprot_mt2_human):
    headers = [("Accept", "application/json"), ("Content-Type", "application/json")]

    res = client.get("/api/v0/uniprots/mt2_human", headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data["description"] == "Metallothionein-2"
    assert data["uniprot_id"] == "MT2_HUMAN"
    assert data["uniprot_acc"] == "P02795"
