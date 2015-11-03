from flask import render_template, request, session
from application import app
import os


@app.route('/<name>')
def content(name):
    return render_template(name, name=name)


@app.route('/')
def index():
    return content('/help.html')
