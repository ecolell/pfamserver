import re


def test_main_app(db, client, main_app_template):
    res = client.get('/')
    assert res.status_code == 200
    data = res.get_data(as_text=True).split('\n')
    diffs = [l for l in zip(main_app_template, data) if l[0] != l[1]]
    assert len(diffs) <= 4
