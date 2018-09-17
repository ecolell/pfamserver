from __future__ import unicode_literals

from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.sql.expression import cast
from sqlalchemy.sql.functions import concat
from sqlalchemy import or_, types
from sqlalchemy.orm import Load
from pfamserver.models import PfamA, PfamARegFullSignificant, Pfamseq, PdbPfamAReg
from pfamserver.extensions import db
from pfamserver.exceptions import SentryIgnoredError
from merry import Merry
from subprocess import Popen as run, PIPE

from pfamserver.services import version_service

merry = Merry()


class PfamServiceError(Exception):
    message = ''

    def __init__(self, message):
        super(PfamServiceError, self).__init__()
        self.message = message


@merry._except(NoResultFound)
def handle_no_result_found(e):
    raise PfamServiceError('PfamA doesn''t exist.')


def get_pfam_acc_from_pfam(code):
    query = db.session.query(PfamA)
    query = query.filter(or_(PfamA.pfamA_acc == code.upper(),
                             PfamA.pfamA_id.ilike(code)))
    return query


@merry._try
def get_pfam(pfam):
    return get_pfam_acc_from_pfam(pfam).one()


def get_sequence_descriptions_from_pfam_with_join_table(pfam, with_pdb):
    subquery = scoped_db.query(PfamA)
    subquery = subquery.filter(or_(PfamA.pfamA_acc == code.upper(),
                                    PfamA.pfamA_id.ilike(code))).distinct().subquery()

    query = scoped_db.query(PfamAPfamseq.pfamseq_id, PfamAPfamseq.pfamA_acc)
    query = query.filter(PfamAPfamseq.pfamA_acc == subquery.c.pfamA_acc)

    if with_pdb:
        query = query.filter(PfamAPfamseq.has_pdb == 1)

    query = query.order_by(PfamAPfamseq.pfamseq_id.asc())
    return query.distinct().all()



def get_sequence_descriptions_from_pfam(pfam, with_pdb):
    subquery = get_pfam_acc_from_pfam(pfam)
    subquery = subquery.distinct().subquery()

    query = db.session.query(concat(Pfamseq.pfamseq_id, '/',
                                    cast(PfamARegFullSignificant.seq_start, types.Unicode), '-',
                                    cast(PfamARegFullSignificant.seq_end, types.Unicode)))
    query = query.join(PfamARegFullSignificant, Pfamseq.pfamseq_acc == PfamARegFullSignificant.pfamseq_acc)
    query = query.filter(PfamARegFullSignificant.pfamA_acc == subquery.c.pfamA_acc)

    if with_pdb:
        subquery2 = db.session.query(PdbPfamAReg)
        subquery2 = subquery2.filter(PdbPfamAReg.pfamA_acc == subquery.c.pfamA_acc).distinct().subquery()
        query = query.filter(PfamARegFullSignificant.pfamseq_acc == subquery2.c.pfamseq_acc)

    query = query.filter(PfamARegFullSignificant.in_full)
    query = query.options(Load(Pfamseq).load_only('pfamseq_id'),
                          Load(PfamARegFullSignificant).load_only("seq_start",
                                                                  "seq_end"))
    query = query.order_by(Pfamseq.pfamseq_id.asc()).distinct()
    results = query.all()
    return [r[0] for r in results]


@merry._try
def get_stockholm_from_pfam(pfam):
    query = get_pfam_acc_from_pfam(pfam)
    query = query.options(Load(PfamA).load_only("pfamA_acc"))
    pfamA_acc = query.one().pfamA_acc
    fetch_call = './esl-afetch'
    cmd = [
        fetch_call,
        "./{version}/Pfam-A.full".format(version=version_service.version()),
        pfamA_acc
    ]
    return run(cmd, stdout=PIPE).communicate()[0]
