from application import app
from flask.ext.restful import Api, Resource
import os
from subprocess import Popen as run, PIPE
from distutils.sysconfig import get_python_lib
from autoupdate import lib_path, db_path
from StringIO import StringIO
from contextlib import closing
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from itertools import chain


api = Api(app)
fetch = '{:s}/hmmer/easel/miniapps/esl-afetch'.format(lib_path)
muscle = '{:s}/muscle/src/muscle -maxiters 1 -diags1 -quiet -sv -distance1 kbit20_3'.format(lib_path)


def get_register(query):
    cmd = [fetch, db_path, query]
    return run(cmd, stdout=PIPE).communicate()[0]


def fill(seqrecord, length):
    seq = seqrecord.seq.__dict__
    seq["_data"] = seq["_data"].ljust(length, '-')
    return seqrecord


def merge(registers):
    pfams = map(lambda reg: StringIO(reg), registers)
    i_msa =  map(lambda pfam: AlignIO.read(pfam, "stockholm"), pfams)
    length = max(map(lambda pfam: pfam.get_alignment_length(), i_msa))
    t_msa = map(lambda msa:
                map(lambda sr: fill(sr, length), msa),
                i_msa)
    seqrecords = list(chain(*t_msa))
    msa = MultipleSeqAlignment(seqrecords)
    map(lambda pfam: pfam.close(), pfams)
    return msa


def realign(msa):
    with closing(StringIO()) as f_tmp:
        count = AlignIO.write(msa, f_tmp, "fasta")
        msa = f_tmp.getvalue()
    msa = run(muscle.split(' '), stdin=PIPE, stdout=PIPE).communicate(input=msa)[0]
    with closing(StringIO()) as f_out:
        with closing(StringIO(msa)) as f_in:
            count = AlignIO.convert(f_in, "fasta", f_out, "stockholm")
        msa = f_out.getvalue() if count else ""
    return msa


def db(query):
    queries = query.split(',')
    registers = map(get_register, queries)
    return registers[0] if len(registers) == 1 else realign(merge(registers))



class QueryAPI(Resource):

    def get(self, query):
        queries = [query, query.upper(), query.capitalize(), query.lower()]

        for q in queries:
            output = db(q)
            if output:
                return {'query': q, 'output': output}
        return {'query': query, 'output': output}


api.add_resource(QueryAPI, '/api/query/<string:query>', endpoint = 'query')
