from __future__ import unicode_literals
from pfamserver.services import sequence_service as service
import pytest


def test_get_pfams_from_sequence(app):
    sequence = 'MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTV' \
               'ERIPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMSMDFQNHLGSCQKCDP' \
               'SCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGCTGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFG' \
               'ATCVKKCPRNYVVTDHGSCVRACGADSYEMEEDGVRKCKKCEGPCRKVCNGIGIGEFKDSLSINATNIKHFKNCTSISGDLHILPVAFRGDSFTHT' \
               'PPLDPQELDILKTVKEITGFLLIQAWPENRTDLHAFENLEIIRGRTKQHGQFSLAVVSLNITSLGLRSLKEISDGDVIISGNKNLCYANTINWKKL' \
               'FGTSGQKTKIISNRGENSCKATGQVCHALCSPEGCWGPEPRDCVSCRNVSRGRECVDKCNLLEGEPREFVENSECIQCHPECLPQAMNITCTGRGP' \
               'DNCIQCAHYIDGPHCVKTCPAGVMGENNTLVWKYADAGHVCHLCHPNCTYGCTGPGLEGCPTNGPKIPSIATGMVGALLLLLVVALGIGLFMRRRH' \
               'IVRKRTLRRLLQERELVEPLTPSGEAPNQALLRILKETEFKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMAS' \
               'VDNPHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAARNVLVKTPQHVKITDFGLAKLLGA' \
               'EEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPK' \
               'FRELIIEFSKMARDPQRYLVIQGDERMHLPSPTDSNFYRALMDEEDMDDVVDADEYLIPQQGFFSSPSTSRTPLLSSLSATSNNSTVACIDRNGLQ' \
               'SCPIKEDSFLQRYSSDPTGALTEDSIDDTFLPVPEYINQSVPKRPAGSVQNPVYHNQPLNPAPSRDPHYQDPHSTAVGNPEYLNTVQPTCVNSTFD' \
               'SPAHWAQKGSHQISLDNPDYQQDFFPKEAKPNGIFKGSTAENAEYLRVAPQSSEFIGA'
    expected_result = [
        {
            "seq_end": 168,
            "pfamA_acc": "PF01030",
            "num_full": 3152,
            "description": "Receptor L domain",
            "seq_start": 57
        },
        {
            "seq_end": 338,
            "pfamA_acc": "PF00757",
            "num_full": 1146,
            "description": "Furin-like cysteine rich region",
            "seq_start": 177
        },
        {
            "seq_end": 481,
            "pfamA_acc": "PF01030",
            "num_full": 3152,
            "description": "Receptor L domain",
            "seq_start": 361
        },
        {
            "seq_end": 637,
            "pfamA_acc": "PF14843",
            "num_full": 1070,
            "description": "Growth factor receptor domain IV",
            "seq_start": 505
        },
        {
            "seq_end": 968,
            "pfamA_acc": "PF07714",
            "num_full": 68047,
            "description": "Protein tyrosine kinase",
            "seq_start": 712
        }
    ]
    with app.app_context():
        assert service.get_pfams_from_sequence(sequence) == expected_result


def test_get_pfams_from_invalid_sequence(app, egfr_human_partial_sequence):
    sequence = egfr_human_partial_sequence
    expected_result = [
        {
            'seq_end': 61,
            'pfamA_acc': 'PF14843',
            'num_full': 1070,
            'description': 'Growth factor receptor domain IV',
            'seq_start': 1
        },
        {
            'seq_end': 392,
            'pfamA_acc': 'PF07714',
            'num_full': 68047,
            'description': 'Protein tyrosine kinase',
            'seq_start': 136
        }
    ]
    with app.app_context():
        assert service.get_pfams_from_sequence(sequence) == expected_result
        sequence = sequence[:30] + '\n' + sequence[31:]
        assert service.get_pfams_from_sequence(sequence) == []


def test_get_pfam_from_pfamacc(db):
    pfamA = service.get_pfam_from_pfamacc('PF00131')
    assert pfamA.description == 'Metallothionein'
    assert pfamA.num_full == 345

    with pytest.raises(service.SequenceServiceError) as exc:
        service.get_pfam_from_pfamacc('invalid_pfam')
        assert exc.value.message == 'PfamA doesn''t exist.'
