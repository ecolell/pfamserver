from application import app
from database import scoped_db as db
from models import classes
from flask_admin import Admin
from flask_admin.contrib.sqla import ModelView


class CompleteView(ModelView):
    column_display_pk = True


admin = Admin(app, name='pfamserver', template_mode='bootstrap3')
for cls in classes:
    admin.add_view(CompleteView(cls, db))
