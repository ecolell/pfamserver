from __future__ import unicode_literals
import json
import pytest


def test_get_protein_sequence(db, client, egfr_human_partial_sequence):
    sequence = egfr_human_partial_sequence
    headers = [('Accept', 'application/json'),
               ('Content-Type', 'application/json')]

    res = client.get('/api/v0/protein_sequences/' + sequence, headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data['query'] == sequence
    results = data['output']
    assert len(results) == 2
    assert results[0]['seq_end'] == 61
    assert results[0]['seq_start'] == 1
    assert results[0]['pfamA_acc'] == 'PF14843'
    assert results[0]['description'] == 'Growth factor receptor domain IV'
    assert 'num_full' in results[0]
    assert results[1]['seq_end'] == 392
    assert results[1]['seq_start'] == 136
    assert results[1]['pfamA_acc'] == 'PF07714'
    assert results[1]['description'] == 'Protein tyrosine kinase'
    assert 'num_full' in results[1]

    sequence = sequence[:30] + '\n' + sequence[31:]
    res = client.get('/api/v0/protein_sequences/' + sequence, headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data['query'] == sequence
    assert data['output'] == []
