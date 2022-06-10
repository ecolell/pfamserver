import pytest
from pfamserver.services import pfam_service as service


def test_get_sequence_descriptions_from_pfam_without_join_table(
    db, pfamseq_pf00131, pfam_a_pfamseq_pf00131, pfamseq_pf01030, pfam_a_pfamseq_pf01030
):
    sequence_descriptions = (
        service.get_sequence_descriptions_from_pfam_without_join_table("pf00131", True)
    )
    assert len(sequence_descriptions) == 6

    sequence_descriptions = (
        service.get_sequence_descriptions_from_pfam_without_join_table("pf00131", False)
    )
    assert len(sequence_descriptions) == 480

    sequence_descriptions = (
        service.get_sequence_descriptions_from_pfam_without_join_table("pf01030", True)
    )
    assert len(sequence_descriptions) == 18

    sequence_descriptions = (
        service.get_sequence_descriptions_from_pfam_without_join_table("pf01030", False)
    )
    assert len(sequence_descriptions) == 4669


def test_get_sequence_descriptions_from_pfam_with_join_table(
    db, pfamseq_pf00131, pfam_a_pfamseq_pf00131, pfamseq_pf01030, pfam_a_pfamseq_pf01030
):
    sequence_descriptions = service.get_sequence_descriptions_from_pfam_with_join_table(
        "pf00131", True
    )
    assert len(sequence_descriptions) == 6

    sequence_descriptions = service.get_sequence_descriptions_from_pfam_with_join_table(
        "pf00131", False
    )
    assert len(sequence_descriptions) == 480

    sequence_descriptions = service.get_sequence_descriptions_from_pfam_with_join_table(
        "pf01030", True
    )
    assert len(sequence_descriptions) == 18

    sequence_descriptions = service.get_sequence_descriptions_from_pfam_with_join_table(
        "pf01030", False
    )
    assert len(sequence_descriptions) == 4669


@pytest.mark.parametrize("table_cache_enabled", [True, False])
def test_get_sequence_descriptions_from_pfam(
    app,
    pfamseq_pf00131,
    pfam_a_pfamseq_pf00131,
    pfamseq_pf01030,
    pfam_a_pfamseq_pf01030,
    table_cache_enabled,
):
    app.config["TABLE_CACHE_ENABLED"] = table_cache_enabled
    sequence_descriptions = service.get_sequence_descriptions_from_pfam("pf00131", True)
    assert len(sequence_descriptions) == 6

    sequence_descriptions = service.get_sequence_descriptions_from_pfam(
        "pf00131", False
    )
    assert len(sequence_descriptions) == 480

    sequence_descriptions = service.get_sequence_descriptions_from_pfam("pf01030", True)
    assert len(sequence_descriptions) == 18

    sequence_descriptions = service.get_sequence_descriptions_from_pfam(
        "pf01030", False
    )
    assert len(sequence_descriptions) == 4669
    app.config["TABLE_CACHE_ENABLED"] = False


def test_get_stockholm_from_pfam_pf00131(
    db, pfam_a_pf00131, mock_pfam_a_pf00131_stockholm
):
    stockholm = service.get_stockholm_from_pfam("pf00131")
    assert len(stockholm) >= 71033  # bytes


def test_get_stockholm_from_pfam_not_mocked(db, pfam_a_pf00131):
    with pytest.raises(service.PfamServiceError) as exc:
        service.get_stockholm_from_pfam("fake_code")
    assert exc.value.message == "PfamA doesn" "t exist."
