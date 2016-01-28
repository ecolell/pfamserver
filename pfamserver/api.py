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
import random
import multiprocessing
from models.version import Version


thread_count = multiprocessing.cpu_count() * 2
print("Working with {:} threads.".format(thread_count))
api = Api(app)
fetch_call = '{:s}/hmmer/easel/miniapps/esl-afetch'.format(lib_path)
muscle_call = '{:s}/muscle/src/muscle -maxiters 1 -diags1 -quiet -sv -distance1 kbit20_3'.format(lib_path)
mafft_call = 'MAFFT_BINARIES={0} {0}/mafft --retree 2 --maxiterate 0 --thread {1} --quiet'.format(lib_path + '/mafft/core', thread_count)


def get_register(query):
    cmd = [fetch_call, db_path, query]
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


def muscle(msa):
    return run(muscle_call.split(' '), stdin=PIPE, stdout=PIPE).communicate(input=msa)[0]


def mafft(msa):
    hash = random.getrandbits(128)
    file_in = '{:}.fasta'.format(hash)
    file_out = '{:}_out.fasta'.format(hash)
    with open(file_in, 'w') as f:
        f.write(msa)
    os.system('{:} {:} > {:}'.format(mafft_call, file_in, file_out))
    with open(file_out, 'r') as f:
        msa = f.read()
    os.system('rm {:} {:}'.format(file_in, file_out))
    return msa


algorithms = {
    "muscle": muscle,
    "mafft": mafft
}


def realign(msa, algorithm):
    with closing(StringIO()) as f_tmp:
        count = AlignIO.write(msa, f_tmp, "fasta")
        msa = f_tmp.getvalue()
    msa = algorithms[algorithm](msa)
    with closing(StringIO()) as f_out:
        with closing(StringIO(msa)) as f_in:
            count = AlignIO.convert(f_in, "fasta", f_out, "stockholm")
        msa = f_out.getvalue() if count else ""
    return msa


def db(query):
    queries = query.split(',')
    registers = map(get_register, queries)
    return registers[0] if len(registers) == 1 else realign(merge(registers), "mafft")



class QueryAPI(Resource):

    def get(self, query):
        queries = [query, query.upper(), query.capitalize(), query.lower()]

        for q in queries:
            output = db(q)
            if output:
                return {'query': q, 'output': output}
        return {'query': query, 'output': output}


api.add_resource(QueryAPI, '/api/query/<string:query>', endpoint = 'query')
