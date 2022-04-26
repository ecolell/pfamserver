from merry import Merry
from pfamserver.extensions import db
from pfamserver.models import PfamA, Uniprot, UniprotRegFull
from sqlalchemy import or_
from sqlalchemy.orm.exc import NoResultFound

merry = Merry()


class UniprotServiceError(Exception):
    message = ""

    def __init__(self, message):
        super().__init__()
        self.message = message


@merry._except(NoResultFound)
def handle_no_result_found(e):
    raise UniprotServiceError("Uniprot doesn" "t exist.")


def query_uniprot(uniprot: str):
    uniprot = uniprot.upper()
    query = db.session.query(Uniprot)
    query = query.join(Uniprot.pfams)
    query = query.filter(
        or_(Uniprot.uniprot_id == uniprot, Uniprot.uniprot_acc == uniprot)
    )
    return query


@merry._try
def get_uniprot(uniprot: str):
    return query_uniprot(uniprot).one()


@merry._try
def get_pfams_from_uniprot(uniprot):
    query = query_uniprot(uniprot)
    query = query.filter(UniprotRegFull.uniprot_acc == Uniprot.uniprot_acc)
    query = query.filter(PfamA.pfamA_acc == UniprotRegFull.pfamA_acc)
    query = query.order_by(UniprotRegFull.seq_start)
    return query.one()
