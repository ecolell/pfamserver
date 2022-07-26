import pytest
from pfamserver.services import pfam_service
from tests.factories import PfamAFactory


@pytest.fixture
def pfam_a_pf01030(db):
    return PfamAFactory(
        pfamA_acc="PF01030",
        pfamA_id="Recep_L_domain",
        num_full=3152,
        description="Receptor L domain",
    )


@pytest.fixture
def pfam_a_pf07714(db):
    return PfamAFactory(
        pfamA_acc="PF07714",
        pfamA_id="Pkinase_Tyr",
        num_full=68047,
        description="Protein tyrosine kinase",
    )


@pytest.fixture
def pfam_a_pf00069(db):
    return PfamAFactory(
        pfamA_acc="PF00069",
        pfamA_id="Pkinase",
        num_full=236455,
        description="Protein kinase domain",
    )


@pytest.fixture
def pfam_a_pf14843(db):
    return PfamAFactory(
        pfamA_acc="PF14843",
        pfamA_id="GF_recep_IV",
        num_full=1070,
        description="Growth factor receptor domain IV",
    )


@pytest.fixture
def pfam_a_pf00757(db):
    return PfamAFactory(
        pfamA_acc="PF00757",
        pfamA_id="Furin-like",
        num_full=1146,
        description="Furin-like cysteine rich region",
    )


@pytest.fixture
def pfam_a_pf00131(db):
    return PfamAFactory(
        pfamA_acc="PF00131",
        pfamA_id="Metallothio",
        num_full=480,
        description="Metallothionein",
    )


@pytest.fixture
def mock_pfam_a_pf00131_stockholm(db, load_input, mocker):
    content = load_input("pf00131.sto")
    m = mocker.patch.object(pfam_service, "get_stockholm_from_pfam")
    m.return_value = content
    return content
