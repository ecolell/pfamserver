from __future__ import unicode_literals

from sqlalchemy.orm.exc import NoResultFound
from flask import current_app
from pfamserver.models import PfamA
from pfamserver.extensions import db
from pfamserver.exceptions import SentryIgnoredError
from merry import Merry
from subprocess import Popen as run, PIPE
import re
import os
from pfamserver.services import version_service
import uuid
from builtins import str as text

merry = Merry()

os.environ["PERL5LIB"] = os.path.abspath('./PfamScan')
os.environ["PATH"] = os.path.abspath('.') + ':' + os.environ["PATH"]


class SequenceServiceError(Exception):
    message = ''

    def __init__(self, message):
        super(SequenceServiceError, self).__init__()
        self.message = message


@merry._except(NoResultFound)
def handle_no_result_found(e):
    raise SequenceServiceError('PfamA doesn''t exist.')


@merry._try
def get_pfam_from_pfamacc(pfam_acc):
    query = db.session.query(PfamA.num_full, PfamA.description)
    query = query.filter(PfamA.pfamA_acc == pfam_acc)
    return query.one()


def is_pfam_match(line):
    return re.match(r"\w+\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+(\w+)\.\d+", line)


def id_generator():
    return text(uuid.uuid4())


def pfamscan(seq):
    hmmdata_path = os.path.abspath(os.path.join('./', version_service.version()))
    pfamscan_bin = os.path.abspath('./PfamScan/pfam_scan.pl')
    tmp_path = os.path.abspath('./tmp')
    pfamscan_call = pfamscan_bin + ' -dir ' + hmmdata_path

    fasta_path = os.path.join(tmp_path, id_generator() + ".fasta")
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


def get_pfams_from_sequence(seq):
    sequence = r'^[AC-IK-Y]*\r*$'
    if re.match(sequence, seq):
        output = parse_pfamscan(pfamscan(seq))
    else:
        output = []
    return output
