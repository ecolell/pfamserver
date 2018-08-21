"""
from application import app, cache
from database import scoped_db
from flask.ext.restless import APIManager
from sqlalchemy import or_, types
from sqlalchemy.orm import Load
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.sql.expression import cast
from sqlalchemy.sql.functions import concat
from models import classes
if classes:
    from models import Uniprot, UniprotRegFull, PfamA, PdbPfamAReg, Pdb, \
        PfamARegFullSignificant, Pfamseq
from flask.ext.restful import Api, Resource
from flask_restful.inputs import boolean
from flask import request
import os
from subprocess import Popen as run, PIPE
from StringIO import StringIO
from contextlib import closing
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from itertools import chain
import random
import multiprocessing
from zlib import compress
from base64 import b64encode
from autoupdate.core import Manager


manager = APIManager(app, flask_sqlalchemy_db=scoped_db)

for cls in classes:
    manager.create_api(cls, methods=['GET'])


thread_count = multiprocessing.cpu_count() * 2
print("Working with {:} threads.".format(thread_count))
api = Api(app)
lib_path = app.config['LIB_PATH']
root_path = Manager().actual_version_path()
fetch_call = '{:s}/hmmer/easel/miniapps/esl-afetch'.format(lib_path)


class StockholmFromPfamAPI(Resource):

    def query(self, query):
        cmd = [fetch_call, root_path + "Pfam-A.full", query]
        return run(cmd, stdout=PIPE).communicate()[0]

    def to_pfam_acc(self, code):
        subquery = scoped_db.query(PfamA)
        subquery = subquery.filter(or_(PfamA.pfamA_acc == code.upper(),
                                       PfamA.pfamA_id.ilike(code)))
        subquery = subquery.options(Load(PfamA).load_only("pfamA_acc"))
        try:
            return subquery.one().pfamA_acc
        except NoResultFound as e:
            return None

    @cache.cached(timeout=3600)
    def get(self, query):
        pfamA_acc = self.to_pfam_acc(query)
        output = ''
        if pfamA_acc:
            output = self.query(pfamA_acc)
        return {'query': pfamA_acc,
                'output': b64encode(compress(output))}




api.add_resource(StockholmFromPfamAPI,
                 '/api/query/stockholm_pfam/<string:query>',
                 endpoint='stockholm_pfam')
"""
