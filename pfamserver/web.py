from flask import render_template, request, session
from application import app, cache
from autoupdate.core import Manager
import os
import json


@cache.cached(timeout=3600)
@app.route('/<name>')
def content(name):
    version = Manager().config['actual_version']
    return render_template(name, name=name, version=version)

@cache.cached(timeout=3600)
@app.route('/')
def index():
    return content('/help.html')
