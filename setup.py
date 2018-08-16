# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import unicode_literals

import os
from pip._internal.req import parse_requirements
from setuptools import setup, find_packages


def read(file_name):
    """ read requirements and create packages list """
    return open(os.path.join(os.path.dirname(__file__), file_name)).read()


# reqs is a list of requirement
# e.g. ['django==1.5.1', 'mezzanine==1.4.6']
reqs = [str(ir.req) for ir in parse_requirements('requirements.txt', session=False)]

setup(
    name='pfamserver',
    version='2.0',
    packages=find_packages(),
    include_package_data=True,
    install_requires=reqs,
    url='pfamserver.leloir.org.ar',
    license='MIT License',
    author='Eloy Adonis Colell',
    author_email='eloy.colell@gmail.com',
    description='PFam Database Server',
    long_description=read('README.md'),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Topic :: Utilities',
        'Programming Language :: Python :: 2.7',
        'License :: MIT License'
    ],
    zip_safe=True,
    platforms=['Linux']
)
