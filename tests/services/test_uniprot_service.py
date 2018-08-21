from __future__ import unicode_literals
from pfamserver.services import uniprot_service as service


def test_get_pfams_from_uniprot(db):
    uniprot = service.get_pfams_from_uniprot('egfr_human')
    assert len(uniprot.pfams) == 6
