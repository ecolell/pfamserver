from __future__ import unicode_literals

from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.sql.expression import cast
from sqlalchemy.sql.functions import concat
from sqlalchemy import or_, types
from sqlalchemy.orm import Load
from flask import current_app
from pfamserver.models import PfamA, PfamARegFullSignificant, Pfamseq, PdbPfamAReg
from pfamserver.extensions import db
from pfamserver.exceptions import SentryIgnoredError
from merry import Merry
from subprocess import Popen as run, PIPE
import re
import os


merry = Merry()


class SequenceServiceError(Exception):
    message = ''

    def __init__(self, message):
        super(SequenceServiceError, self).__init__()
        self.message = message


@merry._except(NoResultFound)
def handle_no_result_found(e):
    raise SequenceServiceError('Sequence doesn''t exist.')


def get_pfam_from_pfamacc(pfam_acc):
    query = db.session.query(PfamA.num_full, PfamA.description)
    query = query.filter(PfamA.pfamA_acc == pfam_acc)
    return query.one()


def is_pfam_match(line):
    return re.match(r"\w+\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+(\w+)\.\d+", line)


def id_generator():
    return '1'


def pfamscan(seq):
    os.environ["PERL5LIB"] = current_app.config['PFAMSCAN_PATH']
    pfamscan_call = ('{:s}/pfam_scan.pl -dir {:s}').format(current_app.config['PFAMSCAN_PATH'], current_app.config['PFAMSCAN_PATH'])

    # fasta_path = "{:s}/PfamScan/P00533.fasta.txt".format(pfamscan_dir)
    fasta_path = os.path.join(current_app.config['TEMP'], id_generator() + ".fasta")
    with open(fasta_path, 'w') as outstream:
        outstream.write(">user_sequence\n" + seq)

    cmd = pfamscan_call.split() + ["-fasta", fasta_path]
    return run(cmd, stdout=PIPE).communicate()[0]


def parse_pfamscan(text):
    matches = [is_pfam_match(line) for line in text.split("\n") if is_pfam_match(line)]
    pfams = [get_pfam_from_pfamacc(m.group(3)) for m in matches]
    return [
        {
            "description": t[1].description,
            "pfamA_acc": t[0].group(3),
            "seq_start": int(t[0].group(1)),
            "seq_end": int(t[0].group(2)),
            "num_full": t[1].num_full
        }
        for t in zip(matches, pfams)]


def get_pfams_from_sequence(sequence):
    seq = sequence.upper().strip()
    sequence = r'^[AC-IK-Y]*\r*$'
    if re.match(sequence, seq):
        output = parse_pfamscan(pfamscan(seq))
    else:
        output = []
    return {
        'query': seq,
        'output': output
    }
