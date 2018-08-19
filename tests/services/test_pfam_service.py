from __future__ import unicode_literals
from pfamserver.services import pfam_service as service


def test_get_sequence_descriptions_from_pfam(db):
    sequence_descriptions = service.get_sequence_descriptions_from_pfam('pf00131', True)
    assert len(sequence_descriptions) == 6

    sequence_descriptions = service.get_sequence_descriptions_from_pfam('pf00131', False)
    assert len(sequence_descriptions) == 345

    sequence_descriptions = service.get_sequence_descriptions_from_pfam('pf01030', True)
    assert len(sequence_descriptions) == 18

    sequence_descriptions = service.get_sequence_descriptions_from_pfam('pf01030', False)
    assert len(sequence_descriptions) == 3152