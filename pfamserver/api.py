from application import app
from flask.ext.restful import Api, Resource
import os
from subprocess import Popen as run, PIPE
from distutils.sysconfig import get_python_lib
from autoupdate import lib_path, db_path


api = Api(app)
fetch = '{:s}/hmmer/easel/miniapps/esl-afetch'.format(lib_path)


def db(query):
    cmd = [fetch, db_path, query]
    return run(cmd, stdout=PIPE).communicate()[0]


class QueryAPI(Resource):

    def get(self, query):
        queries = [query, query.upper(), query.capitalize(), query.lower()]

        for q in queries:
            output = db(q)
            if output:
                return {'query': q, 'output': output}
        return {'query': query, 'output': output}


api.add_resource(QueryAPI, '/api/query/<string:query>', endpoint = 'query')
