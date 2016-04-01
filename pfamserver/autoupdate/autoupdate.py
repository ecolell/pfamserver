# -*- coding: utf-8 -*-
from __future__ import print_function
from distutils.sysconfig import get_python_lib
import json
import socket
import os
import gzip
import shutil
from urllib import urlopen
from contextlib import closing
import subprocess
from progressbar import ProgressBar, ETA, Percentage, RotatingMarker, Bar
import re
from datetime import datetime, timedelta


lib_path = get_python_lib()
#filename = 'uniprot_trembl.dat'
filename = 'uniprot_sprot.dat'
db_path = '{:s}/uniprot_data/{:s}'.format(lib_path, filename)
server = {'path': '/pub/databases/uniprot/current_release/knowledgebase/complete/',
          'url': 'ftp.uniprot.org'}
proto = ['ftp://']
select_proto = None
config_file = '{:s}/uniprot_data/version.json'.format(lib_path)


def save_local_config(config):
    with open(config_file, "w") as f:
        json.dump(config, f)


def get_local_config():
    if not os.path.exists(config_file):
        folder = '/'.join(config_file.split('/')[:-1])
        if not os.path.exists(folder):
            os.mkdir(folder)
        config = {'version': '2010_01',
                  'status': []}
        save_local_config(config)
    with open(config_file, 'r') as f:
        config = json.load(f)
    return config


def silent_remove(filename):
    if os.path.exists(filename):
        os.remove(filename)

@milestone
def clean_old_db(filename):
    silent_remove(filename)
    return True


def get_max_version(protocol):
    url = '{0}{1[url]}{1[path]}reldate.txt'.format(protocol, server)
    with closing(urlopen(url)) as conn:
        if "ftp" in protocol:
            line = conn.read()
            versions = re.findall(r'Release ([0-9\_]*) consists', line)
        elif "http" in protocol:
            # TODO: uniprot doesn't support the http procotol.
            hrefs = filter(lambda l: "href=\"Pfam" in l, conn.readlines())
            versions = map(lambda l: float(re.sub('^.+href="Pfam([0-9\.]+).+$',
                                              r'\1',l)), hrefs)
    return max(versions)


def get_available(config):
    global select_proto
    available = None
    print("->\tUniprot: Get availables versions from {:s}".format(server['url']))
    for p in proto:
        try:
            available = get_max_version(p)
            select_proto = p
            break
        except socket.error:
            print("{:s} failed.".format(p))
            # If the server is offline assume there is no new version.
    if not available:
        print("Discarded new version search.")
        available = config['version']
    return available


def get_versions():
    config = get_local_config()
    remote = get_available(config)
    return config['version'], remote


@milestone
def download(origin, destiny):
    return not subprocess.call("wget -c {} -O {}".format(origin, destiny).split(" "))


@milestone
def download_gziped(remote):
    path_pattern = '{0}{1[url]}{1[path]}{3}.gz'
    origin = path_pattern.format(select_proto, server, remote, filename)
    destiny = '{:s}.gz'.format(db_path)
    clean_old_db(destiny)
    ready = download(origin, destiny)
    return destiny if ready else ready


@milestone
def export(destiny):
    print("->\tUniprot: Extracting database")
    silent_remove(db_path)
    with open(db_path, 'wb') as f_out, gzip.open(destiny, 'rb') as f_in:
        ready = False
        try:
            shutil.copyfileobj(f_in, f_out)
            ready = True
        except Exception:
            pass
    return ready


@milestone
def resume_database():
    print("->\tUniprot: Select the records ID and the Pfam DR (database "
          "references)")
    dt = datetime.now() + timedelta(minutes=3)
    print("  \tThis should be ready for {}".format(dt.time()))
    temp_file = db_path + ".tmp"
    transform = '{} > {}'.format(db_path, temp_file)
    command = 'awk "/^(ID   |DR   Pfam)/{ print $1; }" '
    command += transform
    result = not os.system(command)
    print("->\tUniprot: Updating database")
    if result:
        shutil.move(temp_file, db_path)
    return result


def update():
    local, remote = get_versions()
    if local < remote:
        print("->\tUniprot: Update from {:s} to {:s}".format(local, remote))
        destiny = download_gzied(remote)
        export(destiny)
        resume_database()
        config = {'version': remote, 'status': []}
        save_local_config(config)
    else:
        print("->\tUniprot: The database is updated (version {})".format(local))
