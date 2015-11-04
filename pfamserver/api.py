from application import app
from flask.ext.restful import Api, Resource
import os
from subprocess import Popen as run, PIPE


api = Api(app)


class QueryAPI(Resource):

    def get(self, query):
        cmd = ['./hmmer/binaries/esl-afetch', 'Pfam-A.full', query]
        output = run(cmd, stdout=PIPE).communicate()[0]

        return {'query': query,
                'output': output}


api.add_resource(QueryAPI, '/api/query/<string:query>', endpoint = 'query')
