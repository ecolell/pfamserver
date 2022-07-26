import pytest
from pfamserver.services import uniprot_service as service


def test_get_pfams_from_uniprot(db, uniprot_reg_full_mt2_human):
    uniprot = service.get_pfams_from_uniprot("mt2_human")
    assert len(uniprot.pfams) == 1

    with pytest.raises(service.UniprotServiceError) as exc:
        service.get_pfams_from_uniprot("fake_uniprot")
    assert exc.value.message == "Uniprot doesn" "t exist."


def test_get_uniprot(db, uniprot_reg_full_mt2_human):
    uniprot = service.get_uniprot("mt2_human")
    assert uniprot.description == "Metallothionein-2"
    assert uniprot.uniprot_id == "MT2_HUMAN"
    assert uniprot.uniprot_acc == "P02795"

    with pytest.raises(service.UniprotServiceError) as exc:
        service.get_uniprot("invalid_uniprot")
        assert exc.value.message == "Uniprot doesn" "t exist."
