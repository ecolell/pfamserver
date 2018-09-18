from __future__ import unicode_literals
import json


def test_get_uniprot_pfams(db, client):
    headers = [('Accept', 'application/json'),
               ('Content-Type', 'application/json')]
    res = client.get('/api/v0/uniprots/egfr_human/pfams', headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data['query'] == 'EGFR_HUMAN'
    assert len(data['output']) == 6

    res = client.get('/api/v0/uniprots/invalid_uniprot/pfams', headers=headers)
    assert res.status_code == 400
    data = json.loads(res.get_data(as_text=True))
    assert data['message'] == 'Uniprot doesn''t exist.'


def test_get_reference_sequences_pdbs(db, client):
    headers = [('Accept', 'application/json'),
               ('Content-Type', 'application/json')]
    res = client.get('/api/v0/uniprots/MT2_HUMAN/1-61/pdbs', headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data['query'] == {
        'seq_end': 61,
        'seq_start': 1,
        'uniprot_id': 'MT2_HUMAN'
    }
    assert len(data['output']) == 2


def test_get_uniprots(db, client):
    headers = [('Accept', 'application/json'),
               ('Content-Type', 'application/json')]

    res = client.get('/api/v0/uniprots/egfr_human', headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data['description'] == 'Epidermal growth factor receptor EC=2.7.10.1'
    assert data['uniprot_id'] == 'EGFR_HUMAN'
    assert data['uniprot_acc'] == 'P00533'