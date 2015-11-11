from __future__ import print_function
import urllib2
from distutils.sysconfig import get_python_lib
import json
import socket
import os
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


def silent_remove(filename):
    if os.path.exists(filename):
        os.remove(filename)


def get_max_version():
    conn = urllib2.urlopen('{0[proto]}{0[url]}{0[path]}'.format(ftp))
    versions = map(lambda l: float(l.split(' ')[-1][4:-2]), conn.readlines())
    return max(versions)


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


def save_local_config(config):
    with open(config_file, "w") as f:
        json.dump(config, f)


def get_local_config():
    if not os.path.exists(config_file):
        folder = '/'.join(config_file.split('/')[:-1])
        if not os.path.exists(folder):
            os.mkdir(folder)
        config = {'version': 1.0,
                  'status': []}
        save_local_config(config)
    with open(config_file, 'r') as f:
        config = json.load(f)
    return config


def get_versions(config):
    remote = get_available(config)
    return config['version'], remote


def download_gziped(remote, config):
    origin = '{0[proto]}{0[url]}{0[path]}Pfam{1}/{2}.gz'.format(ftp,
                                                                remote,
                                                                filename)
    destiny = '{:s}.gz'.format(db_path)
    if not 'downloading' in config['status']:
        silent_remove(destiny)
        config['status'].append('downloading')
        save_local_config(config)
    if not 'downloaded' in config['status']:
        download(url=origin, path=destiny)
        config['status'].append('downloaded')
        save_local_config(config)
    return destiny


def export(destiny, config):
    if not 'exported' in config['status']:
        try:
            print("->\tPFam: Extracting database")
            silent_remove(db_path)
            with open(db_path, 'wb') as f_out, gzip.open(destiny, 'rb') as f_in:
                shutil.copyfileobj(f_in, f_out)
            config['status'].append('exported')
            save_local_config(config)
            return True
        except Exception:
            return False


def reindex(config):
    if not 'reindexed' in config['status']:
        print("->\tPFam: Reindexing database")
        from pfamserver.api import fetch
        silent_remove('{}.ssi'.format(db_path))
        print("{0} --index {1}".format(fetch, db_path))
        os.system('{0} --index {1}'.format(fetch, db_path))
        config['status'].append('reindex')
        save_local_config(config)


def update():
    config = get_local_config()
    local, remote = get_versions(config)
    if local < remote:
        print("->\tPFam: Update from {:.1f} to {:.1f}".format(local, remote))
        exported = False
        while not exported:
            destiny = download_gziped(remote, config)
            exported = export(destiny, config)
        reindex(config)
        config = {'version': remote, 'status': []}
        save_local_config(config)
    else:
        print("->\tPFam: The database is updated (version {})".format(local))

