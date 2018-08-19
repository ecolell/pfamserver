from __future__ import unicode_literals

from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm import Load
from pfamserver.models import Uniprot, UniprotRegFull, Pdb, PdbPfamAReg
from pfamserver.extensions import db
from pfamserver.exceptions import SentryIgnoredError
from merry import Merry

merry = Merry()


class SequenceDescriptionServiceError(Exception):
    message = ''

    def __init__(self, message):
        super(SequenceDescriptionServiceError, self).__init__()
        self.message = message


@merry._except(NoResultFound)
def handle_no_result_found(e):
    raise SequenceDescriptionServiceError('PfamA desn''t exists.')


def get_pdbs_from_uniprot_pfam_a_reg(uniprot_id, seq_start, seq_end):
    query = db.session.query(PdbPfamAReg)
    query = query.join(PdbPfamAReg.pdb)
    query = query.join(PdbPfamAReg.uniprot_reg_full)
    query = query.join(UniprotRegFull.uniprot)
    query = query.filter(Uniprot.uniprot_id == uniprot_id,
                         UniprotRegFull.seq_start == seq_start,
                         UniprotRegFull.seq_end == seq_end,
                         UniprotRegFull.uniprot_acc == Uniprot.uniprot_acc,
                         UniprotRegFull.auto_uniprot_reg_full == PdbPfamAReg.auto_uniprot_reg_full,
                         PdbPfamAReg.pdb_id == Pdb.pdb_id)
    query = query.order_by(PdbPfamAReg.pdb_id)
    query = query.order_by(PdbPfamAReg.chain)
    query = query.options(Load(PdbPfamAReg).load_only('pdb_id', 'chain', 'pdb_res_start', 'pdb_res_end', 'pfamA_acc'),
                          Load(Pdb).load_only('title', 'resolution', 'method', 'date', 'author'))
    return query.all()