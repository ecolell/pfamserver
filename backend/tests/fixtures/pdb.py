import pytest
from tests.factories import PdbFactory


@pytest.fixture
def pdb_1mhu(db):
    return PdbFactory(
        pdb_id="1MHU",
        title="""THE THREE-DIMENSIONAL STRUCTURE OF HUMAN [113CD7] METALLOTHIONEIN-2
        IN SOLUTION DETERMINED BY NUCLEAR MAGNETIC RESONANCE SPECTROSCOPY""",
        method="Solution NMR",
        resolution=0,
        author="Braun, W., Messerle, B.A., Schaeffer, A., Vasak, M., Kaegi, J.H.R., Wuthrich, K.",
        date="2011-07-13",
    )


@pytest.fixture
def pdb_2mhu(db):
    return PdbFactory(
        pdb_id="2MHU",
        title="""THE THREE-DIMENSIONAL STRUCTURE OF HUMAN [113CD7] METALLOTHIONEIN-2
        IN SOLUTION DETERMINED BY NUCLEAR MAGNETIC RESONANCE SPECTROSCOPY""",
        method="Solution NMR",
        resolution=0,
        author="Braun, W., Messerle, B.A., Schaeffer, A., Vasak, M., Kaegi, J.H.R., Wuthrich, K.",
        date="2011-07-13",
    )
