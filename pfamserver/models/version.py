from database import db


class Version(db.Model):
    __tablename__ = 'version'

    pfam_release = db.Column(db.Unicode, primary_key=True)
    pfam_release_date = db.Column(db.DateTime)
    swiss_prot_version = db.Column(db.Unicode)
    trembl_version = db.Column(db.Unicode)
    hmmer_version = db.Column(db.Unicode)
    pfamA_coverage = db.Column(db.Float(4,1))
    pfamA_residue_coverage = db.Column(db.Float(4,1))
    number_families = db.Column(db.Integer)
