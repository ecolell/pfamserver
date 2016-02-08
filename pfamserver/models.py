from application import app
from database import engine, scoped_db as db
from sqlalchemy.ext.automap import automap_base
import loader
from sqlalchemy.ext.declarative import declared_attr


class Base(object):

    @declared_attr
    def __tablename__(cls):
        return cls.__name__

    __table_args__ = {'mysql_engine': 'InnoDB'}


Base = automap_base(cls=Base)
Base.query = db.query_property()
Base.prepare(engine, reflect=True)


(
    Version, PfamA, PfamSeq, Uniprot,
    PfamARegSeed, UniprotRegFull, PfamARegFullSignificant, PfamARegFullInsignificant,
    PfamAnnSeq, SecondaryPfamSeqAcc, Evidence,
    MarkupKey, PfamSeqDisulphide, OtherReg, PfamSeqMarkup,
    Architecture, PfamAArchitecture,
    LiteratureReference, GeneOntology, PfamADatabaseLinks, Interpro, PfamALiteratureReference, PfamAInteractions,
    ClanAlignmentAndRelationship, Clan, ClanDatabaseLinks, ClanMembership, ClanLitRef, ClanArchitecture,
    DeadFamily, DeadClan,
    NestedLocations,
    PDB, PDBImage, PDBResidueData, PDBPfamAReg,
    ProteomeRegions, CompleteProteomes, Taxonomy, NCBITaxonomy,
    PfamA2PfamAScoop, PfamA2PfamAHHSearch,
    PfamA_HMM, AlignmentAndTree
) = map(Base.classes.get, loader.tables)
classes = [
    Version, PfamA, PfamSeq, Uniprot,
    PfamARegSeed, UniprotRegFull, PfamARegFullSignificant, PfamARegFullInsignificant,
    PfamAnnSeq, SecondaryPfamSeqAcc, Evidence,
    MarkupKey, PfamSeqDisulphide, OtherReg, PfamSeqMarkup,
    Architecture, PfamAArchitecture,
    LiteratureReference, GeneOntology, PfamADatabaseLinks, Interpro, PfamALiteratureReference, PfamAInteractions,
    ClanAlignmentAndRelationship, Clan, ClanDatabaseLinks, ClanMembership, ClanLitRef, ClanArchitecture,
    DeadFamily, DeadClan,
    NestedLocations,
    PDB, PDBImage, PDBResidueData, PDBPfamAReg,
    ProteomeRegions, CompleteProteomes, Taxonomy, NCBITaxonomy,
    PfamA2PfamAScoop, PfamA2PfamAHHSearch,
    PfamA_HMM, AlignmentAndTree
]
classes = filter(lambda cls: cls, classes)
print "--> Coverage: ", len(classes), "/", len(loader.tables)
