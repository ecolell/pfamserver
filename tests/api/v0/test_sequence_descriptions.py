from __future__ import unicode_literals
import json


def test_get_reference_sequences_pdbs(db, client):
    headers = [('Accept', 'application/json'),
               ('Content-Type', 'application/json')]
    res = client.get('/api/v0/sequence_descriptions/MT2_HUMAN/1-61/pdbs', headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data['query'] == {
        'seq_end': 61,
        'seq_start': 1,
        'uniprot_id': 'MT2_HUMAN'
    }
    assert len(data['output']) == 2
