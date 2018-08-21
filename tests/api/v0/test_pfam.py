from __future__ import unicode_literals
import json


def test_get_pfams_reference_sequences(db, client):
    headers = [('Accept', 'application/json'),
               ('Content-Type', 'application/json')]
    res = client.get('/api/v0/pfams/pf00131/sequence_descriptions', headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data['query'] == 'pf00131'
    assert len(data['output']) == 6

    res = client.get('/api/v0/pfams/pf00131/sequence_descriptions?with_pdb=true', headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data['query'] == 'pf00131'
    assert len(data['output']) == 6

    res = client.get('/api/v0/pfams/pf00131/sequence_descriptions?with_pdb=false', headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data['query'] == 'pf00131'
    assert len(data['output']) == 345

    res = client.get('/api/v0/pfams/PF01030/sequence_descriptions?with_pdb=false', headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data['query'] == 'PF01030'
    assert len(data['output']) == 3152
