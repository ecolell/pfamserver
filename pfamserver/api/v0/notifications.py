# -*- coding: utf-8 -*- #
"""Notifications on API v0."""
from __future__ import unicode_literals
from flask_restplus import Resource

from pfamserver.api.v0 import api, schemas
from pfamserver.services import notification_service


ns = api.namespace('notifications', decorators=[
    api.response(204, "success"),
    api.response(400, "not found")])


@ns.route('')
class ContactUsCreateAPI(Resource):
    """Contact Us Create API."""

    @ns.response(204, "response")
    @ns.doc('Create a new feedback message.')
    @ns.expect(schemas.contact_us_query, validate=True)
    def post(self):
        """Build a contact with our system message."""
        kwargs = schemas.contact_us_query.parse_args()
        notification_service.send_feedback(
            kwargs.get('name', ''),
            kwargs.get('email', ''),
            kwargs.get('phone', ''),
            kwargs.get('message', ''))
        return {}, 204
