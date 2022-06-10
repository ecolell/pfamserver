import json


def test_get_protein_sequence(
    db,
    client,
    egfr_human_partial_sequence,
    mock_pfam_scan_egfr_human,
    uniprot_reg_full_egfr_human,
):
    sequence = egfr_human_partial_sequence
    headers = [("Accept", "application/json"), ("Content-Type", "application/json")]

    res = client.get("/api/v0/protein_sequences/" + sequence, headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data["query"] == sequence
    results = data["output"]
    assert results == [
        {
            "description": "Receptor L domain",
            "pfamA_acc": "PF01030",
            "seq_start": 57,
            "seq_end": 168,
            "num_full": 3152,
        },
        {
            "description": "Furin-like cysteine rich region",
            "pfamA_acc": "PF00757",
            "seq_start": 177,
            "seq_end": 338,
            "num_full": 1146,
        },
        {
            "description": "Receptor L domain",
            "pfamA_acc": "PF01030",
            "seq_start": 361,
            "seq_end": 481,
            "num_full": 3152,
        },
        {
            "description": "Growth factor receptor domain IV",
            "pfamA_acc": "PF14843",
            "seq_start": 505,
            "seq_end": 637,
            "num_full": 1070,
        },
        {
            "description": "Protein tyrosine kinase",
            "pfamA_acc": "PF07714",
            "seq_start": 712,
            "seq_end": 968,
            "num_full": 68047,
        },
    ]

    sequence = sequence[:30] + "\n" + sequence[31:]
    res = client.get("/api/v0/protein_sequences/" + sequence, headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data["query"] == sequence
    assert data["output"] == []
