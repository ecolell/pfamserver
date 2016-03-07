from application import app
from database import scoped_db as db
from models import classes
from flask_admin import Admin
from flask_admin.contrib.sqla import ModelView


specific_attributes = ["pdb_image_sml", "pdb_image", "sequence", "hmm", "logo",
                       "relationship", "alignment", "image_map", "stockholm",
                       "tree", "jtml", "post"]


def view_constructor(cls_in, db):
    attributes = filter(lambda e: e[0][0] != '_', cls_in.__dict__.items())
    blobs = dict(filter(lambda e: "blob" in
                        e[1].type.__class__.__name__.lower(),
                        attributes)).keys()
    related = map(lambda e: e[0], filter(lambda e: e[0][0] != "_" and
                                         e[1].expression.foreign_keys,
                                         attributes))
    not_searchables = related + blobs
    searchables = filter(lambda k: k not in not_searchables,
                         dict(attributes).keys())
    not_filterables = specific_attributes + blobs
    filtereables = filter(lambda k: k not in not_filterables, searchables)
    config = {"column_display_pk": True}
    if blobs:
        config["column_exclude_list"] = tuple(blobs)
    if related:
        config["column_select_related_list"] = tuple(related)
    if searchables:
        config["column_searchable_list"] = tuple(searchables)
    if filtereables:
        config["column_filters"] = tuple(filtereables)
    return type(cls_in.__name__ + "View", (ModelView, ), config)(cls_in, db)


admin = Admin(app, name='pfamserver', template_mode='bootstrap3')
for cls in classes:
    admin.add_view(view_constructor(cls, db))
