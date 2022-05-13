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
    pfam_a_pf00757,
):
    data = [
        [235767053, pfam_a_pf01030, 57, 168, 1],
        [235767055, pfam_a_pf01030, 361, 481, 1],
        [234834757, pfam_a_pf07714, 712, 968, 1],
        [246673367, pfam_a_pf00069, 712, 967, 0],
        [343645098, pfam_a_pf14843, 505, 637, 1],
        [343645098, pfam_a_pf00757, 177, 338, 1],
    ]
    return [
        UniprotRegFullFactory(
            auto_uniprot_reg_full=d[0],
            pfamA=d[1],
            uniprot=uniprot_egfr_human,
            seq_start=d[2],
            seq_end=d[3],
            in_full=d[4],
        )
        for d in data
    ]


@pytest.fixture
def uniprot_reg_full_mt2_human(
    db,
    uniprot_mt2_human,
    pfam_a_pf00131,
):
    return UniprotRegFullFactory(
        auto_uniprot_reg_full=337061979,
        pfamA=pfam_a_pf00131,
        uniprot=uniprot_mt2_human,
        seq_start=1,
        seq_end=61,
        in_full=1,
    )


@pytest.fixture
def uniprot_reg_full_mt3_human(
    db,
    uniprot_mt3_human,
    pfam_a_pf00131,
):
    return UniprotRegFullFactory(
        auto_uniprot_reg_full=337061767,
        pfamA=pfam_a_pf00131,
        uniprot=uniprot_mt3_human,
        seq_start=1,
        seq_end=68,
        in_full=1,
    )
