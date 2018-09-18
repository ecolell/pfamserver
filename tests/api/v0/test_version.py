from __future__ import unicode_literals
import json


def test_get_version(app, client):
    headers = [('Accept', 'application/json'),
               ('Content-Type', 'application/json')]

    res = client.get('/api/v0/version', headers=headers)
    assert res.status_code == 200
    data = json.loads(res.get_data(as_text=True))
    assert data['version'] == 'Pfam31.0'
