import pytest
from tests.factories import UniprotFactory, UniprotRegFullFactory


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
        [344073459, pfam_a_pf14843, 505, 637, 1],
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


@pytest.fixture
def uniprot_reg_full_pfam_a_pf00131(db, pfam_a_pf00131):
    data = [
        (337061767, "P25713", 1, 68, 1),
        (337061779, "P28184", 1, 68, 1),
        (337061979, "P02795", 1, 61, 1),
        (337061993, "P18055", 1, 62, 1),
        (337062085, "P02802", 1, 61, 1),
        (337062186, "P04355", 1, 61, 1),
        (337062415, "P62339", 1, 60, 1),
    ]
    for d in data:
        UniprotFactory(
            uniprot_acc=d[1], uniprot_id=f"{d[1]}_key", description="something"
        )
    return [
        UniprotRegFullFactory(
            auto_uniprot_reg_full=d[0],
            pfamA=pfam_a_pf00131,
            uniprot_acc=d[1],
            seq_start=d[2],
            seq_end=d[3],
            in_full=d[4],
        )
        for d in data
    ]
