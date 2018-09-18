from __future__ import unicode_literals
import json
import pytest


def test_get_protein_sequence(db, client, egfr_human_partial_sequence):
    sequence = egfr_human_partial_sequence
    expected_result = [
        {
            'seq_end': 61,
            'pfamA_acc': 'PF14843',
            'num_full': 1070,
            'description': 'Growth factor receptor domain IV',
            'seq_start': 1
        },
        {
            'seq_end': 392,
            'pfamA_acc': 'PF07714',
            'num_full': 68047,
            'description': 'Protein tyrosine kinase',
            'seq_start': 136
        }
    ]
    headers = [('Accept', 'application/json'),
               ('Content-Type', 'application/json')]

    res = client.get('/api/v0/protein_sequences/' + sequence, headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data['query'] == sequence
    assert data['output'] == expected_result

    sequence = sequence[:30] + '\n' + sequence[31:]
    res = client.get('/api/v0/protein_sequences/' + sequence, headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data['query'] == sequence
    assert data['output'] == []