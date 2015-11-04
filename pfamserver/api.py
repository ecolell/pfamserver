from application import app
from flask.ext.restful import Api, Resource
import os
from subprocess import Popen as run, PIPE


api = Api(app)


def db(query):
    cmd = ['./hmmer/binaries/esl-afetch', 'Pfam-A.full', query]
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
