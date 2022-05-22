import pytest
from pfamserver.services import sequence_service


@pytest.fixture
def egfr_human_complete_sequence():
    return (
        "MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDLSFLKTIQEVAGYVLIALNTV"
        "ERIPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFSNNPALCNVESIQWRDIVSSDFLSNMSMDFQNHLGSCQKCDP"
        "SCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGKSPSDCCHNQCAAGCTGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFG"
        "ATCVKKCPRNYVVTDHGSCVRACGADSYEMEEDGVRKCKKCEGPCRKVCNGIGIGEFKDSLSINATNIKHFKNCTSISGDLHILPVAFRGDSFTHT"
        "PPLDPQELDILKTVKEITGFLLIQAWPENRTDLHAFENLEIIRGRTKQHGQFSLAVVSLNITSLGLRSLKEISDGDVIISGNKNLCYANTINWKKL"
        "FGTSGQKTKIISNRGENSCKATGQVCHALCSPEGCWGPEPRDCVSCRNVSRGRECVDKCNLLEGEPREFVENSECIQCHPECLPQAMNITCTGRGP"
        "DNCIQCAHYIDGPHCVKTCPAGVMGENNTLVWKYADAGHVCHLCHPNCTYGCTGPGLEGCPTNGPKIPSIATGMVGALLLLLVVALGIGLFMRRRH"
        "IVRKRTLRRLLQERELVEPLTPSGEAPNQALLRILKETEFKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMAS"
        "VDNPHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAARNVLVKTPQHVKITDFGLAKLLGA"
        "EEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPK"
        "FRELIIEFSKMARDPQRYLVIQGDERMHLPSPTDSNFYRALMDEEDMDDVVDADEYLIPQQGFFSSPSTSRTPLLSSLSATSNNSTVACIDRNGLQ"
        "SCPIKEDSFLQRYSSDPTGALTEDSIDDTFLPVPEYINQSVPKRPAGSVQNPVYHNQPLNPAPSRDPHYQDPHSTAVGNPEYLNTVQPTCVNSTFD"
        "SPAHWAQKGSHQISLDNPDYQQDFFPKEAKPNGIFKGSTAENAEYLRVAPQSSEFIGA"
    )


@pytest.fixture
def egfr_human_partial_sequence():
    return (
        "DNCIQCAHYIDGPHCVKTCPAGVMGENNTLVWKYADAGHVCHLCHPNCTYGCTGPGLEGCPTNGPKIPSIATGMVGALLLLLVVALGIGLFMRRRH"
        "IVRKRTLRRLLQERELVEPLTPSGEAPNQALLRILKETEFKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMAS"
        "VDNPHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAARNVLVKTPQHVKITDFGLAKLLGA"
        "EEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPK"
        "FRELIIEFSKMARDPQRYLVIQGDERMHLPSPTDSNFYRALMDEEDMDDVVDADEYLIPQQGFFSSPSTSRTPLLSSLSATSNNSTVACIDRNGLQ"
        "SCPIKEDSFLQRYSSDPTGALTEDSIDDTFLPVPEYINQSVPKRPAGSVQNPVYHNQPLNPAPSRDPHYQDPHSTAVGNPEYLNTVQPTCVNSTFD"
        "SPAHWAQKGSHQISLDNPDYQQDFFPKEAKPNGIFKGSTAENAEYLRVAPQSSEFIGA"
    )


@pytest.fixture
def shared_pfam_scan_header():
    return (
        b"# pfam_scan.pl,  run at Sun May 22 05:00:30 2022\n#\n# Copyright (c) 2009 Genome Research Ltd\n"
        b"# Freely distributed under the GNU \n# General Public License\n#\n# Authors: Jaina Mistry (jain"
        b"a@ebi.ac.uk), \n#          Rob Finn (rdf@ebi.ac.uk)\n#\n# This is free software; you can redist"
        b"ribute it and/or modify it under\n# the terms of the GNU General Public License as published by"
        b" the Free Software\n# Foundation; either version 2 of the License, or (at your option) any late"
        b"r version.\n# This program is distributed in the hope that it will be useful, but WITHOUT\n# AN"
        b"Y WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS\n# FOR A PARTICULAR"
        b" PURPOSE. See the GNU General Public License for more\n# details.\n#\n# You should have receive"
        b"d a copy of the GNU General Public License along with\n# this program. If not, see <http://www."
        b"gnu.org/licenses/>. \n# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ="
        b" = = = =\n#      query sequence file: /home/pfamserver/stage/tmp/9f6c44c6-10fd-4bcf-a415-a135a8"
        b"c9793b.fasta\n#        searching against: /home/pfamserver/stage/Pfam35.0/Pfam-A.hmm, with cut "
        b"off --cut_ga\n#    resolve clan overlaps: on\n#     predict active sites: off\n# = = = = = = = "
        b"= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n#\n# <seq id> <alignment sta"
        b"rt> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hm"
        b"m end> <hmm length> <bit score> <E-value> <significance> <clan>\n"
    )


@pytest.fixture
def mock_pfam_scan_egfr_human(mocker, shared_pfam_scan_header):
    a = (
        b"\nuser_sequence     57    167 "
        b"    57    168 PF01030.27  Recep_L_domain    Repeat     1   109   110    107.3   4.6e-31   1 CL0"
        b"022   \nuser_sequence    185    338    177    338 PF00757.23  Furin-like        Domain    10   "
        b"149   149    101.3   4.9e-29   1 CL0547   \nuser_sequence    361    480    361    481 PF01030.2"
        b"7  Recep_L_domain    Repeat     1   109   110     95.2   2.7e-27   1 CL0022   \nuser_sequence  "
        b"  505    636    505    637 PF14843.9   GF_recep_IV       Domain     1   131   132    157.2   1."
        b"9e-46   1 CL0547   \nuser_sequence    713    965    712    968 PF07714.20  PK_Tyr_Ser-Thr    Do"
        b"main     2   256   259    294.2   8.1e-88   1 CL0016   \n"
    )
    m = mocker.patch.object(sequence_service, "pfamscan")
    m.return_value = shared_pfam_scan_header + a


@pytest.fixture
def mock_pfam_scan_partial_egfr_human(mocker, shared_pfam_scan_header):
    a = (
        b"\nuser_sequence      1     60      1     61 PF14843.9   GF_recep_IV       Domain    73   131   13"
        b"2     65.6   4.1e-18   1 CL0547   \nuser_sequence    137    389    136    392 PF07714.20  PK_Tyr_"
        b"Ser-Thr    Domain     2   256   259    295.9   2.4e-88   1 CL0016   \n"
    )
    m = mocker.patch.object(sequence_service, "pfamscan")
    m.return_value = shared_pfam_scan_header + a
