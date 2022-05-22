import os
import re
import uuid
from builtins import str as text
from subprocess import PIPE  # nosec
from subprocess import Popen as run  # nosec

from merry import Merry
from pfamserver.extensions import db
from pfamserver.models import PfamA
from pfamserver.services import version_service
from sqlalchemy.orm.exc import NoResultFound

merry = Merry()

os.environ["PERL5LIB"] = (
    os.path.abspath("./Pfam35.0/PfamScan")
    + ":"
    + os.path.abspath("./Pfam35.0")
    + ":"
    + os.path.abspath(".")
)  # TODO
os.environ["PATH"] = (
    os.path.abspath("./Pfam35.0")
    + ":"
    + os.path.abspath(".")
    + ":"
    + os.environ["PATH"]
)


class SequenceServiceError(Exception):
    message = ""

    def __init__(self, message):
        super().__init__()
        self.message = message


@merry._except(NoResultFound)
def handle_no_result_found(e):
    raise SequenceServiceError("PfamA doesn" "t exist.")


@merry._try
def get_pfam_from_pfamacc(pfam_acc):
    query = db.session.query(PfamA.num_full, PfamA.description)
    query = query.filter(PfamA.pfamA_acc == pfam_acc)
    return query.one()


def is_pfam_match(line):
    return re.match(r"\w+\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+(\w+)\.\d+", line)


def id_generator():
    return text(uuid.uuid4())


def guarantee_tmp_folder(tmp):
    if not os.path.exists(tmp):
        os.makedirs(tmp)


def pfamscan(seq):
    HMMDATA_PATH = os.path.abspath(os.path.join("./", version_service.version()))
    PFAMSCAN_BIN = os.path.join(HMMDATA_PATH, "PfamScan", "pfam_scan.pl")
    TMP_PATH = os.path.abspath("./tmp")
    PFAMSCAN_BASE_CALL = PFAMSCAN_BIN + " -dir " + HMMDATA_PATH
    guarantee_tmp_folder(TMP_PATH)
    fasta_path = os.path.join(TMP_PATH, id_generator() + ".fasta")
    with open(fasta_path, "w") as outstream:
        outstream.write(">user_sequence\n" + seq)

    cmd = PFAMSCAN_BASE_CALL.split() + ["-fasta", fasta_path]
    return run(cmd, stdout=PIPE).communicate()[0]  # nosec


def parse_pfamscan(text):
    print(text)
    lines = text.decode("unicode_escape").split("\n")
    matches = [is_pfam_match(line) for line in lines if is_pfam_match(line)]
    pfams = [get_pfam_from_pfamacc(m.group(3)) for m in matches]
    return [
        {
            "description": t[1].description,
            "pfamA_acc": t[0].group(3),
            "seq_start": int(t[0].group(1)),
            "seq_end": int(t[0].group(2)),
            "num_full": t[1].num_full,
        }
        for t in zip(matches, pfams)
    ]


def get_pfams_from_sequence(seq):
    sequence = r"^[AC-IK-Y]*\r*$"
    if re.match(sequence, seq):
        output = parse_pfamscan(pfamscan(seq))
    else:
        output = []
    return output
