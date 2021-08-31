# -*- coding: utf-8 -*- #
"""Notifications Service."""
from __future__ import unicode_literals
from flask_mail import Message
from pfamserver.extensions import mail
from pfamserver.exceptions import SentryIgnoredError
from flask import current_app


class NotificationServiceError(SentryIgnoredError):
    """Class to notify a notification error. It is ignored by Sentry."""

    message = ''

    def __init__(self, message):
        """
        Initialize the exception with the error message.

            :param self: the exception instance.
            :param message: a service message.
            :type message: str
        """
        super(NotificationServiceError, self).__init__()
        self.message = message


def send_feedback(name, email, phone, message):
    """
    Send an email to our feedback email, with the user name, a message, a phone, and an email to respond.

        :param name: the user name.
        :param email: the user email.
        :param phone: the user phone.
        :param message: the user message.
    """
    if not name or not email or not phone or not message:
        raise NotificationServiceError('You forgot to provide a field of the form.')
    subject = 'A new PFAMSERVER message from {name} ({phone})'.format(
        name=name,
        phone=phone)
    msg = Message(recipients=[current_app.config.get('FEEDBACK_EMAIL')],
                  reply_to=email,
                  body=message,
                  subject=subject)
    mail.send(msg)
