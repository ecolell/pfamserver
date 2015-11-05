from __future__ import print_function
from ftplib import FTP
from distutils.sysconfig import get_python_lib
import json
import socket
import os
from subprocess import call
from homura import download
import gzip
import shutil


lib_path = get_python_lib()
db_path = '{:s}/pfam_data/Pfam-A.full'.format(lib_path)
ftp = {'proto': 'ftp://',
       'path': '/pub/databases/Pfam/releases/',
       'url': 'ftp.ebi.ac.uk'}
filename = 'Pfam-A.full'
config_file = '{:s}/pfam_data/PfamA_version.json'.format(lib_path)


def get_max_version():
    conn = FTP(ftp['url'])
    conn.login()
    path = ftp['path']
    versions = conn.nlst(path)
    version = max(map(lambda f: float(f[len(path)+4:]), versions))
    conn.close()
    return version


def get_available(config):
    available = None
    print("->\tPFam: Get availables versions from {:s}".format(ftp['url']))
    while not available:
        try:
            available = get_max_version()
        except socket.error:
            print("Failed")
            # If the server is offline assume there is no new version.
            available = config['version']
    return available


def save_local_version(config):
    with open(config_file, "w") as f:
        json.dump(config, f)


def get_local_version():
    if not os.path.exists(config_file):
        folder = '/'.join(config_file.split('/')[:-1])
        if not os.path.exists(folder):
            os.mkdir(folder)
        config = {'version': 1.0}
        save_local_version(config)
    with open(config_file, 'r') as f:
        config = json.load(f)
    return config


def get_versions(config):
    remote = get_available(config)
    return config['version'], remote


def download_gziped(remote):
    origin = '{0[proto]}{0[url]}{0[path]}Pfam{1}/{2}.gz'.format(ftp,
                                                                remote,
                                                                filename)
    destiny = '{:s}.gz'.format(db_path)
    download(url=origin, path=destiny)
    return destiny


def export(destiny):
    with open(db_path, 'rb') as f_in, gzip.open(destiny, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)


def reindex():
    from pfamserver.api import fetch
    os.remove('{}.ssi'.format(db_path))
    call('{0} --index {1}'.format(fecth, db_path))


def update():
    config = get_local_version()
    local, remote = get_versions(config)
    if local < remote:
        print("->\tPFam: Update from {:.1f} to {:.1f}".format(local, remote))
        destiny = download_gziped(remote)
        export(destiny)
        reindex()
        config['version'] = remote
        save_local_version(config)
