from application import app
from flask.ext.restful import Api, Resource
import os
from subprocess import Popen as run, PIPE
from distutils.sysconfig import get_python_lib


api = Api(app)
fetch = '{:s}/hmmer/easel/miniapps/esl-afetch'.format(get_python_lib())


def db(query):
    cmd = [fetch, 'Pfam-A.full', query]
    return run(cmd, stdout=PIPE).communicate()[0]


class QueryAPI(Resource):

    def get(self, query):
        queries = [query, query.capitalize(), query.upper(), query.lower()]

        for q in queries:
            output = db(q)
            if output:
                return {'query': q, 'output': output}
        return {'query': query, 'output': output}


api.add_resource(QueryAPI, '/api/query/<string:query>', endpoint = 'query')
