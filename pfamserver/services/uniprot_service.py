from __future__ import unicode_literals

from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy import or_

from pfamserver.models import Uniprot, UniprotRegFull, PfamA
from pfamserver.extensions import db, cache
from pfamserver.exceptions import SentryIgnoredError
from merry import Merry

merry = Merry()


class UniprotServiceError(Exception):
    message = ''

    def __init__(self, message):
        super(UniprotServiceError, self).__init__()
        self.message = message


@merry._except(NoResultFound)
def handle_no_result_found(e):
    raise UniprotServiceError('Job already Fixed.')


def get_pfams_from_uniprot(uniprot):
    uniprot = uniprot.upper()
    query = db.session.query(Uniprot)
    query = query.join(Uniprot.pfams)
    query = query.filter(or_(
        Uniprot.uniprot_id == uniprot,
        Uniprot.uniprot_acc == uniprot))
    query = query.filter(UniprotRegFull.uniprot_acc == Uniprot.uniprot_acc)
    query = query.filter(PfamA.pfamA_acc == UniprotRegFull.pfamA_acc)
    query = query.order_by(UniprotRegFull.seq_start)
    return query.first()
