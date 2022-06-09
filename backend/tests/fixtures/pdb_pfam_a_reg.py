import pytest
from tests.factories import PdbFactory, PdbPfamARegFactory


@pytest.fixture
def pdb_pfam_a_reg_pf00131(db, pfam_a_pf00131, uniprot_reg_full_pfam_a_pf00131):
    # auto_uniprot, pdb_id, pfamseq_acc, chain, res_start, res_end
    data = [
        [337061767, "2FJ4", "P25713", "A", 32, 68],
        [337061767, "2F5H", "P25713", "A", 32, 68],
        [337061767, "2FJ5", "P25713", "A", 32, 68],
        [337061779, "1JI9", "P28184", "A", 32, 68],
        [337061979, "1MHU", "P02795", "A", 31, 61],
        [337061979, "2MHU", "P02795", "A", 1, 30],
        [337061993, "2MRB", "P18055", "A", 1, 30],
        [337061993, "1MRB", "P18055", "A", 31, 61],
        [337062085, "1DFT", "P02802", "A", 1, 30],
        [337062085, "1DFS", "P02802", "A", 1, 31],
        [337062186, "2MRT", "P04355", "A", 1, 30],
        [337062186, "4MT2", "P04355", "A", 1, 61],
        [337062186, "1MRT", "P04355", "A", 31, 61],
        [337062415, "1M0J", "P62339", "A", 8, 35],
        [337062415, "1M0G", "P62339", "A", 37, 66],
    ]
    not_repeated_pdb = {d[1] for d in data}
    for d in not_repeated_pdb:
        PdbFactory(pdb_id=d)
    yield [
        PdbPfamARegFactory(
            auto_uniprot_reg_full=d[0],
            pdb_id=d[1],
            chain=d[3],
            pfamA_acc=pfam_a_pf00131.pfamA_acc,
            pfamseq_acc=d[2],
            pdb_res_start=d[4],
            pdb_res_end=d[5],
        )
        for d in data
    ]


@pytest.fixture
def pdb_pfam_a_reg_pf01030(db, pfam_a_pf01030, uniprot_reg_full_pfam_a_pf01030):
    # auto_uniprot, pdb_id, pfamseq_acc, chain, res_start, res_end
    data = [
        (235766120, "5CUS", "P21860", "B", 36, 148),
        (235766120, "4P59", "P21860", "A", 55, 167),
        (235766120, "3P11", "P21860", "A", 36, 148),
        (235766120, "4LEO", "P21860", "C", 36, 148),
        (235766120, "1M6B", "P21860", "B", 36, 148),
        (235766120, "1M6B", "P21860", "A", 36, 148),
        (235766120, "5CUS", "P21860", "A", 36, 148),
        (235766120, "5CUS", "P21860", "D", 38, 148),
        (235766120, "5CUS", "P21860", "C", 36, 148),
        (235766123, "5CUS", "P21860", "B", 334, 455),
        (235766123, "4P59", "P21860", "A", 353, 474),
        (235766123, "3P11", "P21860", "A", 334, 455),
        (235766123, "4LEO", "P21860", "C", 334, 455),
        (235766123, "1M6B", "P21860", "B", 334, 455),
        (235766123, "1M6B", "P21860", "A", 334, 455),
        (235766123, "5CUS", "P21860", "A", 334, 455),
        (235766123, "5CUS", "P21860", "D", 334, 455),
        (235766123, "5CUS", "P21860", "C", 334, 455),
        (235767053, "5WB7", "P00533", "A", 33, 144),
        (235767053, "3NJP", "P00533", "A", 33, 144),
        (235767053, "1NQL", "P00533", "A", 33, 144),
        (235767053, "1IVO", "P00533", "A", 33, 144),
        (235767053, "3QWQ", "P00533", "A", 33, 144),
        (235767053, "4UIP", "P00533", "A", 33, 144),
        (235767053, "1MOX", "P00533", "B", 33, 144),
        (235767053, "5WB8", "P00533", "D", 33, 144),
        (235767053, "3NJP", "P00533", "B", 33, 144),
        (235767053, "5WB7", "P00533", "B", 33, 144),
        (235767053, "1IVO", "P00533", "B", 33, 144),
        (235767053, "5XWD", "P00533", "A", 33, 144),
        (235767053, "4UV7", "P00533", "A", 33, 144),
        (235767053, "4KRO", "P00533", "A", 33, 144),
        (235767053, "5WB7", "P00533", "C", 33, 144),
        (235767053, "5WB8", "P00533", "A", 33, 144),
        (235767053, "4KRP", "P00533", "A", 33, 132),
        (235767053, "1MOX", "P00533", "A", 33, 144),
        (235767053, "5WB7", "P00533", "D", 33, 144),
        (235767053, "1YY9", "P00533", "A", 33, 144),
        (235767055, "5WB7", "P00533", "A", 337, 457),
        (235767055, "6B3S", "P00533", "B", 337, 457),
        (235767055, "4KRM", "P00533", "G", 337, 457),
        (235767055, "4KRL", "P00533", "A", 337, 457),
        (235767055, "3NJP", "P00533", "A", 337, 457),
        (235767055, "1IVO", "P00533", "A", 337, 457),
        (235767055, "1NQL", "P00533", "A", 337, 457),
        (235767055, "5SX5", "P00533", "M", 337, 457),
        (235767055, "3QWQ", "P00533", "A", 337, 457),
        (235767055, "4UIP", "P00533", "A", 337, 457),
        (235767055, "6B3S", "P00533", "E", 337, 457),
        (235767055, "4KRM", "P00533", "C", 337, 457),
        (235767055, "1MOX", "P00533", "B", 337, 457),
        (235767055, "5SX4", "P00533", "N", 337, 457),
        (235767055, "3B2U", "P00533", "S", 337, 457),
        (235767055, "5WB8", "P00533", "D", 337, 457),
        (235767055, "3NJP", "P00533", "B", 337, 457),
        (235767055, "6B3S", "P00533", "I", 337, 457),
        (235767055, "4KRM", "P00533", "E", 337, 457),
        (235767055, "3B2U", "P00533", "M", 337, 457),
        (235767055, "4KRM", "P00533", "A", 337, 457),
        (235767055, "6B3S", "P00533", "A", 337, 457),
        (235767055, "5WB7", "P00533", "B", 337, 457),
        (235767055, "1IVO", "P00533", "B", 337, 457),
        (235767055, "3B2U", "P00533", "P", 337, 457),
        (235767055, "3C09", "P00533", "A", 337, 457),
        (235767055, "3B2U", "P00533", "E", 337, 457),
        (235767055, "5XWD", "P00533", "A", 337, 457),
        (235767055, "3B2U", "P00533", "I", 337, 457),
        (235767055, "5SX5", "P00533", "N", 337, 457),
        (235767055, "3C09", "P00533", "D", 337, 457),
        (235767055, "4UV7", "P00533", "A", 337, 457),
        (235767055, "4KRO", "P00533", "A", 337, 457),
        (235767055, "4KRM", "P00533", "I", 337, 457),
        (235767055, "3P0Y", "P00533", "A", 337, 457),
        (235767055, "3B2V", "P00533", "A", 337, 457),
        (235767055, "3B2U", "P00533", "B", 337, 457),
        (235767055, "5WB7", "P00533", "C", 337, 457),
        (235767055, "5WB8", "P00533", "A", 337, 457),
        (235767055, "5SX4", "P00533", "M", 337, 457),
        (235767055, "4KRP", "P00533", "A", 337, 457),
        (235767055, "4KRM", "P00533", "K", 337, 457),
        (235767055, "1MOX", "P00533", "A", 337, 457),
        (235767055, "5WB7", "P00533", "D", 337, 457),
        (235767055, "1YY9", "P00533", "A", 337, 457),
        (235767055, "3B2U", "P00533", "V", 337, 457),
        (235767055, "3B2U", "P00533", "A", 337, 457),
        (235767604, "3U9U", "Q15303", "E", 30, 142),
        (235767604, "3U7U", "Q15303", "B", 30, 142),
        (235767604, "3U9U", "Q15303", "F", 30, 142),
        (235767604, "3U7U", "Q15303", "A", 30, 142),
        (235767604, "3U7U", "Q15303", "E", 30, 142),
        (235767604, "2AHX", "Q15303", "B", 30, 142),
        (235767604, "2AHX", "Q15303", "A", 30, 142),
        (235767604, "3U7U", "Q15303", "D", 30, 142),
        (235767604, "3U2P", "Q15303", "A", 30, 142),
        (235767604, "3U7U", "Q15303", "C", 30, 142),
        (235767604, "3U7U", "Q15303", "F", 30, 142),
        (235767607, "3U9U", "Q15303", "E", 333, 453),
        (235767607, "3U7U", "Q15303", "B", 333, 453),
        (235767607, "3U9U", "Q15303", "F", 333, 453),
        (235767607, "3U7U", "Q15303", "A", 333, 453),
        (235767607, "3U7U", "Q15303", "E", 333, 453),
        (235767607, "2AHX", "Q15303", "B", 333, 453),
        (235767607, "2AHX", "Q15303", "A", 333, 453),
        (235767607, "3U7U", "Q15303", "D", 333, 453),
        (235767607, "3U2P", "Q15303", "A", 333, 453),
        (235767607, "3U7U", "Q15303", "C", 333, 453),
        (235767607, "3U7U", "Q15303", "F", 333, 453),
        (235770615, "5U8Q", "P08069", "A", 21, 131),
        (235770615, "1IGR", "P08069", "A", 21, 131),
        (235770615, "5U8R", "P08069", "A", 21, 131),
        (235770626, "5U8Q", "P08069", "A", 322, 436),
        (235770626, "1IGR", "P08069", "A", 322, 436),
        (235770626, "5U8R", "P08069", "A", 322, 436),
        (235771374, "3LTG", "P04412", "A", 29, 140),
        (235771374, "3I2T", "P04412", "A", 29, 140),
        (235771374, "3LTF", "P04412", "C", 29, 140),
        (235771374, "3LTF", "P04412", "A", 29, 140),
        (235771374, "3LTG", "P04412", "C", 29, 140),
        (235771377, "3LTG", "P04412", "A", 320, 449),
        (235771377, "3I2T", "P04412", "A", 320, 449),
        (235771377, "3LTF", "P04412", "C", 320, 449),
        (235771377, "3LTF", "P04412", "A", 320, 449),
        (235771377, "3LTG", "P04412", "C", 320, 449),
        (235772598, "6CEB", "P06213", "A", 25, 137),
        (235772598, "2HR7", "P06213", "B", 25, 137),
        (235772598, "3W11", "P06213", "E", 25, 137),
        (235772598, "6CE9", "P06213", "A", 25, 137),
        (235772598, "5KQV", "P06213", "F", 74, 137),
        (235772598, "4XSS", "P06213", "E", 25, 137),
        (235772598, "6CE7", "P06213", "A", 25, 137),
        (235772598, "4XST", "P06213", "E", 25, 137),
        (235772598, "6CEB", "P06213", "B", 25, 137),
        (235772598, "2HR7", "P06213", "A", 25, 137),
        (235772598, "3W12", "P06213", "E", 25, 137),
        (235772598, "3W13", "P06213", "E", 25, 137),
        (235772598, "5J3H", "P06213", "E", 25, 137),
        (235772598, "5KQV", "P06213", "E", 74, 137),
        (235772598, "6CE7", "P06213", "B", 25, 137),
        (235772598, "6CE9", "P06213", "B", 25, 137),
        (235772598, "4OGA", "P06213", "E", 25, 137),
        (235772598, "4ZXB", "P06213", "E", 25, 137),
        (235772605, "6CEB", "P06213", "A", 332, 446),
        (235772605, "2HR7", "P06213", "B", 332, 446),
        (235772605, "6CE9", "P06213", "A", 332, 446),
        (235772605, "5KQV", "P06213", "F", 25, 446),
        (235772605, "6CE7", "P06213", "A", 332, 446),
        (235772605, "6CEB", "P06213", "B", 332, 446),
        (235772605, "2HR7", "P06213", "A", 332, 446),
        (235772605, "5KQV", "P06213", "E", 25, 446),
        (235772605, "6CE7", "P06213", "B", 332, 446),
        (235772605, "6CE9", "P06213", "B", 332, 446),
        (235772605, "4ZXB", "P06213", "E", 332, 446),
        (235773899, "1LK2", "P15208", "P", 1, 8),
        (235774544, "2A91", "P04626", "A", 31, 152),
        (235774544, "3MZW", "P04626", "A", 30, 151),
        (235774544, "3N85", "P04626", "A", 30, 151),
        (235774544, "1S78", "P04626", "A", 30, 151),
        (235774544, "3WSQ", "P04626", "A", 30, 151),
        (235774544, "4HRM", "P04626", "A", 30, 151),
        (235774544, "3BE1", "P04626", "A", 30, 151),
        (235774544, "3WLW", "P04626", "A", 30, 151),
        (235774544, "3H3B", "P04626", "B", 30, 151),
        (235774544, "3H3B", "P04626", "A", 30, 151),
        (235774544, "5K33", "P04626", "C", 30, 151),
        (235774544, "1S78", "P04626", "B", 30, 151),
        (235774544, "4HRM", "P04626", "C", 30, 151),
        (235774544, "5KWG", "P04626", "C", 30, 151),
        (235774544, "3WLW", "P04626", "B", 30, 151),
        (235774544, "1N8Z", "P04626", "C", 30, 151),
        (235774544, "4HRL", "P04626", "C", 30, 151),
        (235774544, "5MY6", "P04626", "A", 52, 173),
        (235774544, "6ATT", "P04626", "A", 30, 151),
        (235774548, "2A91", "P04626", "A", 345, 465),
        (235774548, "3MZW", "P04626", "A", 344, 464),
        (235774548, "3N85", "P04626", "A", 344, 464),
        (235774548, "1S78", "P04626", "A", 344, 464),
        (235774548, "3WSQ", "P04626", "A", 344, 464),
        (235774548, "3BE1", "P04626", "A", 344, 464),
        (235774548, "3WLW", "P04626", "A", 344, 464),
        (235774548, "5K33", "P04626", "C", 344, 464),
        (235774548, "1S78", "P04626", "B", 344, 464),
        (235774548, "5KWG", "P04626", "C", 344, 464),
        (235774548, "3WLW", "P04626", "B", 344, 464),
        (235774548, "1N8Z", "P04626", "C", 344, 464),
        (235774548, "5MY6", "P04626", "A", 366, 486),
        (235774548, "6ATT", "P04626", "A", 344, 464),
        (235775188, "1N8Y", "P06494", "C", 30, 152),
        (235775191, "1N8Y", "P06494", "C", 345, 465),
    ]
    not_repeated_pdb = {d[1] for d in data}
    for d in not_repeated_pdb:
        PdbFactory(pdb_id=d)
    yield [
        PdbPfamARegFactory(
            auto_uniprot_reg_full=d[0],
            pdb_id=d[1],
            chain=d[3],
            pfamA_acc=pfam_a_pf01030.pfamA_acc,
            pfamseq_acc=d[2],
            pdb_res_start=d[4],
            pdb_res_end=d[5],
        )
        for d in data
    ]
