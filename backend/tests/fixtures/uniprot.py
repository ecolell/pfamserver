import pytest
from tests.factories import UniprotFactory


@pytest.fixture
def uniprot_egfr_human():
    return UniprotFactory(
        uniprot_acc="P00533",
        uniprot_id="EGFR_HUMAN",
        description="Epidermal growth factor receptor EC=2.7.10.1",
    )


@pytest.fixture
def uniprot_mt2_human():
    return UniprotFactory(
        uniprot_acc="P02795", uniprot_id="MT2_HUMAN", description="Metallothionein-2"
    )
