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

merry = Merry()


class PfamServiceError(Exception):
    message = ''

    def __init__(self, message):
        super(PfamServiceError, self).__init__()
        self.message = message


@merry._except(NoResultFound)
def handle_no_result_found(e):
    raise PfamServiceError('PfamA desn''t exists.')


def get_sequence_descriptions_from_pfam(pfam, with_pdb):
    # icode = "%{:}%".format(code)
    subquery = db.session.query(PfamA)
    subquery = subquery.filter(or_(PfamA.pfamA_acc == pfam.upper(),
                                   PfamA.pfamA_id.ilike(pfam))).distinct().subquery()

    # query = db.session.query(UniprotRegFull, Uniprot, PdbPfamAReg)
    # query = query.filter(UniprotRegFull.pfamA_acc == subquery.c.pfamA_acc)
    # query = query.filter(UniprotRegFull.auto_uniprot_reg_full == PdbPfamAReg.auto_uniprot_reg_full)
    # query = query.filter(UniprotRegFull.uniprot_acc == Uniprot.uniprot_acc)

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