import pytest


@pytest.fixture(scope="function")
def main_app_template():
    return """<!DOCTYPE html>
<html lang="en">

<head>
    <link rel="shortcut icon" href="{{ url_for('static', filename='favicon.ico') }}">
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta name="description" content="">
    <meta name="author" content="">

    <title>pfamserver</title>

</head>

<body id="page-top">
<div>
    <h2>Examples</h2>
    <div>
        <div>
            <a href="/api/v0/uniprots/egfr_human">1. Uniprot description</a>
        </div>
        <div>
            <a href="/api/v0/uniprots/egfr_human/pfams">2. Available pfam for a uniprot</a>
        </div>
        <div>
            <a href="/api/v0/pfams/pf00131">3. Pfam description</a>
        </div>
        <div>
            <a href="/api/v0/pfams/pf00131/sequence_descriptions?with_pdb=false">4. Pfam known sequences</a>
        </div>
        <div>
            <a href="/api/v0/pfams/pf00131/sequence_descriptions">5. Pfam known sequences with known PDB</a>
        </div>
        <div>
            <a href="/api/v0/pfams/pf00131/stockholm">6. Pfam MSA in stockholm format</a>
        </div>
        <div>
            <a href="/api/v0/uniprots/egfr_human/505-637/pdbs">7. vailable pdbs for a sequence</a>
        </div>
        <div>
            <a href="/api/v0/protein_sequences/DNCIQCAHYIDGPHCVKTCPAGVMGENNTLVWKYADAGHVCHLCHPNCTYGCTGPGLEGCPTNGPKIPSIATGMVGALLLLLVVALGIGLFMRRRHIVRKRTLRRLLQERELVEPLTPSGEAPNQALLRILKETEFKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVDNPHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAARNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPKFRELIIEFSKMARDPQRYLVIQGDERMHLPSPTDSNFYRALMDEEDMDDVVDADEYLIPQQGFFSSPSTSRTPLLSSLSATSNNSTVACIDRNGLQSCPIKEDSFLQRYSSDPTGALTEDSIDDTFLPVPEYINQSVPKRPAGSVQNPVYHNQPLNPAPSRDPHYQDPHSTAVGNPEYLNTVQPTCVNSTFDSPAHWAQKGSHQISLDNPDYQQDFFPKEAKPNGIFKGSTAENAEYLRVAPQSSEFIGA">8. Pfam from a partial protein sequence</a>
        </div>
        <div>
            <a href="/api/v0/version">9. Version</a>
        </div>
    </div>
</div>
<div>
    <h2>Mantainers</h2>
    <div>
        <div>
            <a href="https://www.linkedin.com/in/eloy-adonis-colell-43902720/?locale=en_US">Eloy Colell</a>
        </div>
        <div>
            <a href="mailto:cmb@leloir.org.ar">Cristina Marino Buslje</a>
        </div>
        <div>
            <a href="mailto:fsimonetti@leloir.org.ar">Franco Simonetti</a>
        </div>
        <div>
            <a href="mailto:jiserte@leloir.org.ar">Javier Iserte</a>
        </div>
    </div>
</div>
</body>
</html>
""".split(
        "\n"
    )
