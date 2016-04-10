from flask import render_template, request, session
from application import app
from autoupdate.core import Manager
import os
import json


@app.route('/<name>')
def content(name):
    version = Manager().config['actual_version']
    return render_template(name, name=name, version=version)

@app.route('/')
def index():
    return content('/help.html')
