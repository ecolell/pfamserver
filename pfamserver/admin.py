from application import app
from database import db
from models.version import Version
from models.pfamA import PfamA
from flask_admin import Admin
from flask_admin.contrib.sqla import ModelView


class CompleteView(ModelView):
    column_display_pk = True

admin = Admin(app, name='pfamserver', template_mode='bootstrap3')
admin.add_view(CompleteView(Version, db.session))
admin.add_view(CompleteView(PfamA, db.session))
