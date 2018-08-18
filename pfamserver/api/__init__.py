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


class SequenceDescriptionFromPfamAPI(Resource):

    def get_descriptions(self, code, with_pdb):
        #icode = "%{:}%".format(code)
        subquery = scoped_db.query(PfamA)
        subquery = subquery.filter(or_(PfamA.pfamA_acc == code.upper(),
                                       PfamA.pfamA_id.ilike(code))).distinct().subquery()

        #query = scoped_db.query(UniprotRegFull, Uniprot, PdbPfamAReg)
        #query = query.filter(UniprotRegFull.pfamA_acc == subquery.c.pfamA_acc)
        #query = query.filter(UniprotRegFull.auto_uniprot_reg_full == PdbPfamAReg.auto_uniprot_reg_full)
        #query = query.filter(UniprotRegFull.uniprot_acc == Uniprot.uniprot_acc)

        query = scoped_db.query(concat(Pfamseq.pfamseq_id, '/',
                                       cast(PfamARegFullSignificant.seq_start, types.Unicode), '-',
                                       cast(PfamARegFullSignificant.seq_end, types.Unicode)))
        query = query.join(PfamARegFullSignificant, Pfamseq.pfamseq_acc == PfamARegFullSignificant.pfamseq_acc)
        query = query.filter(PfamARegFullSignificant.pfamA_acc == subquery.c.pfamA_acc)

        if with_pdb:
            subquery2 = scoped_db.query(PdbPfamAReg)
            subquery2 = subquery2.filter(PdbPfamAReg.pfamA_acc == subquery.c.pfamA_acc).distinct().subquery()
            query = query.filter(PfamARegFullSignificant.pfamseq_acc == subquery2.c.pfamseq_acc)

        query = query.filter(PfamARegFullSignificant.in_full)
        query = query.options(Load(Pfamseq).load_only('pfamseq_id'),
                              Load(PfamARegFullSignificant).load_only("seq_start",
                                                                      "seq_end"))
        query = query.order_by(Pfamseq.pfamseq_id.asc())
        return query.distinct().all()

    @cache.memoize(timeout=3600)
    def get(self, query):
        with_pdb = boolean(request.args.get('with_pdb', 'true'))
        response = {'query': query, 'with_pdb': with_pdb}
        output = self.get_descriptions(query, with_pdb)
        if output:
            response['output'] = [o[0] for o in output]
            response['size'] = len(response['output'])
        return response


class PdbFromSequenceDescriptionAPI(Resource):

    def query(self, uniprot_id, seq_start, seq_end):
        query = scoped_db.query(Uniprot, UniprotRegFull, PdbPfamAReg, Pdb)
        query = query.filter(Uniprot.uniprot_id == uniprot_id,
                             UniprotRegFull.seq_start == seq_start,
                             UniprotRegFull.seq_end == seq_end,
                             UniprotRegFull.uniprot_acc == Uniprot.uniprot_acc,
                             UniprotRegFull.auto_uniprot_reg_full == PdbPfamAReg.auto_uniprot_reg_full,
                             PdbPfamAReg.pdb_id == Pdb.pdb_id)
        query = query.order_by(PdbPfamAReg.pdb_id)
        query = query.order_by(PdbPfamAReg.chain)
        query = query.options(Load(PdbPfamAReg).load_only("pdb_id", "chain", "pdb_res_start", "pdb_res_end"),
                              Load(UniprotRegFull).load_only("pfamA_acc"),
                              Load(Pdb).load_only("title", "resolution", "method", "date", "author"))
        return query.all()

    def serialize(self, element):
        authors = element.Pdb.author.split(',')
        return {
            'pdb_id': element.PdbPfamAReg.pdb_id,
            'chain': element.PdbPfamAReg.chain,
            'pdb_res_start': element.PdbPfamAReg.pdb_res_start,
            'pdb_res_end': element.PdbPfamAReg.pdb_res_end,
            'pfamA_acc': element.UniprotRegFull.pfamA_acc,
            'title': element.Pdb.title,
            'resolution': float(element.Pdb.resolution),
            'method': element.Pdb.method,
            'author': (authors[0] + ' et. al.'
                       if len(authors) > 2 else ''),
            'date': element.Pdb.date
        }

    @cache.cached(timeout=3600)
    def get(self, query):
        uniprot_id, seq_start, seq_end = query.split(',')
        response = {
            'query': {
                'uniprot_id': uniprot_id,
                'seq_start': seq_start,
                'seq_end': seq_end}}
        output = self.query(uniprot_id, seq_start, seq_end)
        if output:
            response['output'] = map(self.serialize, output)
        return response


api.add_resource(StockholmFromPfamAPI,
                 '/api/query/stockholm_pfam/<string:query>',
                 endpoint='stockholm_pfam')
api.add_resource(SequenceDescriptionFromPfamAPI,
                 '/api/query/sequencedescription_pfam/<string:query>',
                 endpoint='sequencedescription_pfam')
api.add_resource(PdbFromSequenceDescriptionAPI,
                 '/api/query/pdb_sequencedescription/<string:query>',
                 endpoint='pdb_sequencedescription')
"""