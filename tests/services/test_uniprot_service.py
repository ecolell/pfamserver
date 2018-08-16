from pfamserver.services import uniprot_service as service


def test_get_pfams_from_uniprot(db):
    pfams = service.get_pfams_from_uniprot('egfr_human')
    assert len(pfams) == 6