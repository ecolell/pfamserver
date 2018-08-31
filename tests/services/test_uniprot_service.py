from __future__ import unicode_literals
from pfamserver.services import uniprot_service as service
import pytest


def test_get_pfams_from_uniprot(db):
    uniprot = service.get_pfams_from_uniprot('egfr_human')
    assert len(uniprot.pfams) == 6

    with pytest.raises(service.UniprotServiceError) as exc:
        service.get_pfams_from_uniprot('fake_uniprot')
    assert exc.value.message == 'Uniprot doesn''t exist.'
