from __future__ import unicode_literals
import pytest


@pytest.fixture
def current_version(app):
    return "Pfam32.0"