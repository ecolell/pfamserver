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
        [pfam_a_pf01030, 57, 168, 1],
        [pfam_a_pf01030, 361, 481, 1],
        [pfam_a_pf07714, 712, 968, 1],
        [pfam_a_pf00069, 712, 967, 0],
        [pfam_a_pf14843, 505, 637, 1],
        [pfam_a_pf00757, 177, 338, 1],
    ]
    return [
        UniprotRegFullFactory(
            auto_uniprot_reg_full=1,
            pfamA=d[0],
            uniprot=uniprot_egfr_human,
            seq_start=d[1],
            seq_end=d[2],
            in_full=d[3],
        )
        for d in data
    ]


@pytest.fixture
def uniprot_reg_full_mt2_human_section1(
    db,
    uniprot_mt2_human,
    pfam_a_pf00131,
):
    return UniprotRegFullFactory(
        auto_uniprot_reg_full=2,
        pfamA=pfam_a_pf00131,
        uniprot=uniprot_mt2_human,
        seq_start=57,
        seq_end=168,
        in_full=1,
    )


@pytest.fixture
def uniprot_reg_full_mt2_human_section2(
    db,
    uniprot_mt2_human,
    pfam_a_pf00131,
):
    return UniprotRegFullFactory(
        auto_uniprot_reg_full=2,
        pfamA=pfam_a_pf00131,
        uniprot=uniprot_mt2_human,
        seq_start=361,
        seq_end=481,
        in_full=1,
    )
