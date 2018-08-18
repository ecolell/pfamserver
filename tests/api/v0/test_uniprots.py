from __future__ import unicode_literals
import json


def test_job_request_create_post(db, client):
    headers = [('Accept', 'application/json'),
               ('Content-Type', 'application/json')]
    res = client.get('/api/v0/uniprots/egfr_human/pfams', headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data['query'] == 'EGFR_HUMAN'
    assert len(data['output']) == 6