def test_main_app(client, main_app_template):
    res = client.get("/")
    assert res.status_code == 200
    data = res.get_data(as_text=True).split("\n")
    diffs = [line for line in zip(main_app_template, data) if line[0] != line[1]]
    assert len(diffs) <= 4
