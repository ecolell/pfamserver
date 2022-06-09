import json

import pytest


def test_get_pfams(db, client, pfam_a_pf00131, pfam_a_pf01030):
    headers = [("Accept", "application/json"), ("Content-Type", "application/json")]

    res = client.get("/api/v0/pfams/pf00131", headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data["pfamA_acc"] == "PF00131"
    assert data["pfamA_id"] == "Metallothio"
    assert data["description"] == "Metallothionein"
    assert data["num_full"] >= 345

    res = client.get("/api/v0/pfams/PF01030", headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data["pfamA_acc"] == "PF01030"
    assert data["pfamA_id"] == "Recep_L_domain"
    assert data["description"] == "Receptor L domain"
    assert data["num_full"] >= 3152


@pytest.mark.parametrize("table_cache_enabled", [True, False])
def test_get_pfams_pf00131_reference_sequences(
    app, client, pfamseq_pf00131, pfam_a_pfamseq_pf00131, table_cache_enabled
):
    headers = [("Accept", "application/json"), ("Content-Type", "application/json")]

    app.config["TABLE_CACHE_ENABLED"] = table_cache_enabled
    res = client.get("/api/v0/pfams/pf00131/sequence_descriptions", headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data["query"] == "pf00131"
    assert len(data["output"]) == 6

    res = client.get(
        "/api/v0/pfams/pf00131/sequence_descriptions?with_pdb=true", headers=headers
    )
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data["query"] == "pf00131"
    assert len(data["output"]) == 6

    res = client.get(
        "/api/v0/pfams/pf00131/sequence_descriptions?with_pdb=false", headers=headers
    )
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data["query"] == "pf00131"
    assert len(data["output"]) == 480

    app.config["TABLE_CACHE_ENABLED"] = False


@pytest.mark.parametrize("table_cache_enabled", [True, False])
def test_get_pfams_pf01030_reference_sequences(
    app, client, pfamseq_pf01030, pfam_a_pfamseq_pf01030, table_cache_enabled
):
    headers = [("Accept", "application/json"), ("Content-Type", "application/json")]

    app.config["TABLE_CACHE_ENABLED"] = table_cache_enabled
    res = client.get(
        "/api/v0/pfams/PF01030/sequence_descriptions?with_pdb=false", headers=headers
    )
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data["query"] == "PF01030"
    assert len(data["output"]) == 4669
    app.config["TABLE_CACHE_ENABLED"] = False


def test_get_pfams_stockholm(db, client, pfam_a_pf00131, mock_pfam_a_pf00131_stockholm):
    headers = [("Accept", "application/json"), ("Content-Type", "application/json")]
    res = client.get("/api/v0/pfams/pf00131/stockholm", headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data["query"] == "pf00131"
    assert len(data["output"]) >= 14896


def test_get_pfams_stockholm_invalid(db, client):
    headers = [("Accept", "application/json"), ("Content-Type", "application/json")]
    res = client.get("/api/v0/pfams/invalid/stockholm", headers=headers)
    assert res.status_code == 400
    data = json.loads(res.get_data(as_text=True))
    assert data["message"] == "PfamA doesn" "t exist."
