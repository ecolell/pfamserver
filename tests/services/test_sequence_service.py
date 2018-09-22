from __future__ import unicode_literals
from pfamserver.services import sequence_service as service
import pytest


def test_get_pfams_from_sequence(app, egfr_human_complete_sequence):
    with app.app_context():
        results = service.get_pfams_from_sequence(egfr_human_complete_sequence)
        assert len(results) == 5

        assert results[0]['seq_end'] == 168
        assert results[0]['seq_start'] == 57
        assert results[0]['pfamA_acc'] == 'PF01030'
        assert results[0]['description'] == 'Receptor L domain'
        assert 'num_full' in results[0]

        assert results[1]['seq_end'] == 338
        assert results[1]['seq_start'] == 177
        assert results[1]['pfamA_acc'] == 'PF00757'
        assert results[1]['description'] == 'Furin-like cysteine rich region'
        assert 'num_full' in results[1]

        assert results[2]['seq_end'] == 481
        assert results[2]['seq_start'] == 361
        assert results[2]['pfamA_acc'] == 'PF01030'
        assert results[2]['description'] == 'Receptor L domain'
        assert 'num_full' in results[2]

        assert results[3]['seq_end'] == 637
        assert results[3]['seq_start'] == 505
        assert results[3]['pfamA_acc'] == 'PF14843'
        assert results[3]['description'] == 'Growth factor receptor domain IV'
        assert 'num_full' in results[3]

        assert results[4]['seq_end'] == 968
        assert results[4]['seq_start'] == 712
        assert results[4]['pfamA_acc'] == 'PF07714'
        assert results[4]['description'] == 'Protein tyrosine kinase'
        assert 'num_full' in results[4]


def test_get_pfams_from_invalid_sequence(app, egfr_human_partial_sequence):
    sequence = egfr_human_partial_sequence
    with app.app_context():
        results = service.get_pfams_from_sequence(sequence)
        assert len(results) == 2
        assert results[0]['seq_end'] == 61
        assert results[0]['seq_start'] == 1
        assert results[0]['pfamA_acc'] == 'PF14843'
        assert results[0]['description'] == 'Growth factor receptor domain IV'
        assert 'num_full' in results[0]
        assert results[1]['seq_end'] == 392
        assert results[1]['seq_start'] == 136
        assert results[1]['pfamA_acc'] == 'PF07714'
        assert results[1]['description'] == 'Protein tyrosine kinase'
        assert 'num_full' in results[1]

        sequence = sequence[:30] + '\n' + sequence[31:]
        assert service.get_pfams_from_sequence(sequence) == []


def test_get_pfam_from_pfamacc(db):
    pfamA = service.get_pfam_from_pfamacc('PF00131')
    assert pfamA.description == 'Metallothionein'
    assert pfamA.num_full >= 345

    with pytest.raises(service.SequenceServiceError) as exc:
        service.get_pfam_from_pfamacc('invalid_pfam')
        assert exc.value.message == 'PfamA doesn''t exist.'
