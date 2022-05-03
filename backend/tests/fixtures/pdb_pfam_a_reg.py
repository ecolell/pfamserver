import pytest
from tests.factories import PdbPfamARegFactory


@pytest.fixture
def pdb_pfam_a_reg_mt2_human_1mhu(
    db, pfam_a_pf00131, uniprot_reg_full_mt2_human_section1, pdb_1mhu
):
    PdbPfamARegFactory(
        auto_uniprot_reg_full=uniprot_reg_full_mt2_human_section1.auto_uniprot_reg_full,
        uniprot_reg_full=uniprot_reg_full_mt2_human_section1,
        pdb=pdb_1mhu,
        chain="A",
        pfamA_acc=pfam_a_pf00131.pfamA_acc,
        pfamseq_acc="123",
        pdb_res_start=31,
        pdb_res_end=61,
    )


@pytest.fixture
def pdb_pfam_a_reg_mt2_human_2mhu(
    db, pfam_a_pf00131, uniprot_reg_full_mt2_human_section2, pdb_2mhu
):
    PdbPfamARegFactory(
        auto_uniprot_reg_full=uniprot_reg_full_mt2_human_section2.auto_uniprot_reg_full,
        uniprot_reg_full=uniprot_reg_full_mt2_human_section2,
        pdb=pdb_2mhu,
        chain="A",
        pfamA_acc=pfam_a_pf00131.pfamA_acc,
        pfamseq_acc="123",
        pdb_res_start=1,
        pdb_res_end=30,
    )
