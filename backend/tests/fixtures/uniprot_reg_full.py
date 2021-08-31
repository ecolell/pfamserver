from __future__ import unicode_literals
import pytest
from tests.factories import UniprotRegFullFactory


@pytest.fixture
def uniprot_reg_full_egfr_human(
    db,
    uniprot_egfr_human,
    pfam_a_pf01030,
    pfam_a_pf07714,
    pfam_a_pf00069,
    pfam_a_pf14843,
    pfam_a_pf00757
):
    data = [
        [pfam_a_pf01030, 57, 168, 1],
        [pfam_a_pf01030, 361, 481, 1],
        [pfam_a_pf07714, 712, 968, 1],
        [pfam_a_pf00069, 712, 967, 0],
        [pfam_a_pf14843, 505, 637, 1],
        [pfam_a_pf00757, 177, 338, 1]
    ]
    for d in data:
        UniprotRegFullFactory(
            pfamA=d[0],
            uniprot=uniprot_egfr_human,
            seq_start=d[1],
            seq_end=d[2],
            in_full=d[3]
        )
    return uniprot_egfr_human
