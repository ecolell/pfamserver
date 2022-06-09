import pytest
from tests.factories import PdbFactory, PdbPfamARegFactory


@pytest.fixture
def pdb_pfam_a_reg_pf00131(db, pfam_a_pf00131, uniprot_reg_full_pfam_a_pf00131):
    # auto_uniprot, pdb_id, pfamA_Acc, pfamseq_acc, chain, res_start, res_end
    data = [
        [337061767, "2FJ4", "PF00131", "P25713", "A", 32, 68],
        [337061767, "2F5H", "PF00131", "P25713", "A", 32, 68],
        [337061767, "2FJ5", "PF00131", "P25713", "A", 32, 68],
        [337061779, "1JI9", "PF00131", "P28184", "A", 32, 68],
        [337061979, "1MHU", "PF00131", "P02795", "A", 31, 61],
        [337061979, "2MHU", "PF00131", "P02795", "A", 1, 30],
        [337061993, "2MRB", "PF00131", "P18055", "A", 1, 30],
        [337061993, "1MRB", "PF00131", "P18055", "A", 31, 61],
        [337062085, "1DFT", "PF00131", "P02802", "A", 1, 30],
        [337062085, "1DFS", "PF00131", "P02802", "A", 1, 31],
        [337062186, "2MRT", "PF00131", "P04355", "A", 1, 30],
        [337062186, "4MT2", "PF00131", "P04355", "A", 1, 61],
        [337062186, "1MRT", "PF00131", "P04355", "A", 31, 61],
        [337062415, "1M0J", "PF00131", "P62339", "A", 8, 35],
        [337062415, "1M0G", "PF00131", "P62339", "A", 37, 66],
    ]
    for d in data:
        PdbFactory(pdb_id=d[1])
    yield [
        PdbPfamARegFactory(
            auto_uniprot_reg_full=d[0],
            pdb_id=d[1],
            chain=d[4],
            pfamA_acc=d[2],
            pfamseq_acc=d[3],
            pdb_res_start=d[5],
            pdb_res_end=d[6],
        )
        for d in data
    ]
