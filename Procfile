backend: flask run
backend_wsgi: uwsgi --socket 0.0.0.0:8000 --protocol=http -w wsgi --processes=2 --threads=2
webpack: npm run start
web: gunicorn manage:app
release: flask db upgrade && flask data algorithms load 
check_safety: pip freeze | safety check --stdin
