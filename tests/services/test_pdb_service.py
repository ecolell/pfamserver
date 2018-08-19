from __future__ import unicode_literals
from pfamserver.services import pdb_service as service


def test_get_pdbs_from_uniprot_pfam_a_reg(db):
    sequence_descriptions = service.get_pdbs_from_uniprot_pfam_a_reg('MT2_HUMAN', 1, 61)
    assert len(sequence_descriptions) == 2

    sequence_descriptions = service.get_pdbs_from_uniprot_pfam_a_reg('MT3_HUMAN', 1, 68)
    assert len(sequence_descriptions) == 3
