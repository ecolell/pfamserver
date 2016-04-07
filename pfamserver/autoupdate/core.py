import json
from application import app
import socket
import gzip
import shutil
from urllib import urlopen
from contextlib import closing
import subprocess
from progressbar import ProgressBar, ETA, Percentage, RotatingMarker, Bar
import re
import os


proto = ['ftp://', 'http://']
select_proto = None


class Manager(object):
    root_path = app.config['ROOT_PATH']
    config_file = '{:s}version.json'.format(root_path)

    def __init__(self):
        self.server = {'path': '/pub/databases/Pfam/releases/',
                       'url': 'ftp.ebi.ac.uk'}

    @property
    def config(self):
        if not os.path.exists(Manager.config_file):
            folder = '/'.join(Manager.config_file.split('/')[:-1])
            if not os.path.exists(folder):
                os.mkdir(folder)
            config = {
                'versions': {
                    '1.0': {
                        'status': []
                    }
                },
                'actual_version': '1.0'}
            self.config = config
        with open(Manager.config_file, 'r') as f:
            config = json.load(f)
        return config

    @config.setter
    def config(self, config):
        with open(Manager.config_file, "w") as f:
            json.dump(config, f)

    @staticmethod
    def milestone(function):
        extension = ''
        if isinstance(function, tuple):
            function, extension = function
        manager = Manager()
        def wrapper(*args):
            self = args[0]
            config = manager.config
            version = (config["versions"][self.version]
                      if self.version in config["versions"] else {"status": []})
            params = ''
            if len(args) > 1:
                params = '-' + '_'.join(args[1:])
            name = function.__name__ + extension + params
            pending = name not in version['status']
            ready = False
            if pending:
                try:
                    ready = function(*args)
                except Exception, e:
                    print "*---->", e
                    return False
                if ready:
                    version['status'].append(name)
                    config["versions"][self.version] = version
                    manager.config = config
            return not pending or ready
        return wrapper

    def version_path(self, version):
        return "{:s}{:s}/".format(Manager.root_path, version.version)

    def get_max_version(self, protocol):
        url = '{0}{1[url]}{1[path]}'.format(protocol, self.server)
        with closing(urlopen(url)) as conn:
            if "ftp" in protocol:
                versions = map(lambda l: float(l.split(' ')[-1][4:-2]),
                               conn.readlines())
            elif "http" in protocol:
                hrefs = filter(lambda l: "href=\"Pfam" in l, conn.readlines())
                versions = map(lambda l: float(re.sub('^.+href="Pfam([0-9\.]+).+$',
                                                      r'\1',l)), hrefs)
        return max(versions)

    def get_available(self, config):
        global select_proto
        available = None
        print("->\tPFam: Get availables versions from {:s}".format(self.server['url']))
        for p in proto:
            try:
                available = self.get_max_version(p)
                select_proto = p
                break
            except socket.error:
                print("{:s} failed.".format(p))
                # If the server is offline assume there is no new version.
        if not available:
            print("Discarded new version search.")
            available = config['actual_version']
        return available

    def actual_version_path(self):
        return self.version_path(Version(self.config['actual_version'], self))

    def get_versions(self):
        config = self.config
        remote = self.get_available(config)
        return config['actual_version'], remote

    def update(self):
        local, remote = self.get_versions()
        if float(local) < remote:
            remote = "{:.1f}".format(remote)
            old = Version(local, self)
            new = Version(remote, self)
            print("->\tPFam: Update from {:} to {:}".format(local, remote))
            new.prepare()
            # new.upgrade()
            old.remove()
        else:
            print("->\tPFam: The database is updated (version {})".format(local))

    def silent_remove(self, filename):
        if os.path.exists(filename):
            os.remove(filename)



class Version(object):

    def __init__(self, version, manager):
        self.version = version
        self.manager = manager
        self.calculate_paths()

    def calculate_paths(self):
        self.path = self.manager.version_path(self)
        self.downloads = ["Pfam-A.full.gz", "database_files/"]
        self.downloaded = dict(map(lambda k: (k, self.path + k),
                                   self.downloads))
        self.extracted = {
            "Pfam-A.full.gz": self.path + "Pfam-A.full"
        }

    def download(self, filename):
        if not os.path.exists(self.path):
            os.mkdir(self.path)
        folder = "P" if "/" in filename else "O"
        path_pattern = '{0}{1[url]}{1[path]}Pfam{2}/{3}'
        origin = path_pattern.format(select_proto, self.manager.server,
                                     self.version, filename)
        destiny = self.downloaded[filename]
        command = "wget -c -r {} -{} {}".format(origin, folder, destiny)
        return not subprocess.call(command.split(" "))

    @Manager.milestone
    def download_pfama_full(self):
        return self.download('Pfam-A.full.gz')

    @Manager.milestone
    def export(self):
        print("->\tPFam: Extracting database")
        filename = "Pfam-A.full.gz"
        with open(self.extracted[filename], 'wb') as f_out:
            with gzip.open(self.downloaded[filename], 'rb') as f_in:
                ready = False
                try:
                    shutil.copyfileobj(f_in, f_out)
                    ready = True
                except Exception, e:
                    pass
        return ready

    @Manager.milestone
    def clean_accession_numbers(self):
        print("->\tPFam: Spliting Pfam IDs and Pfam revision")
        filename = "Pfam-A.full.gz"
        temp_file = self.extracted[filename] + ".tmp"
        fsize = lambda f: os.stat(f).st_size
        max_value = fsize(self.extracted[filename]) * 1.10
        with open(temp_file, 'wb') as f_out:
            with open(self.extracted[filename],"r") as f_in:
                widgets = [Percentage(), ' ', Bar(marker='#'), ' ', ETA()]
                with ProgressBar(widgets=widgets,
                                 max_value=max_value) as progress:
                    for line in f_in:
                        new_line = re.sub(r'^(#=GF AC   [A-Z0-9]+)\.(.+)$',
                                          r'\1\n#=GF DC   Revision: \2', line)
                        f_out.write(new_line)
                        progress.update(fsize(temp_file))
            print("->\tPFam: Updating Pfam IDs")
            shutil.move(temp_file, self.extracted[filename])
        return True

    @Manager.milestone
    def reindex(self):
        print("->\tPFam: Reindexing database")
        from pfamserver.api import fetch_call
        filename = "Pfam-A.full.gz"
        index_file = '{}.ssi'.format(self.extracted[filename])
        self.manager.silent_remove(index_file)
        print("{0} --index {1}".format(fetch_call, self.extracted[filename]))
        os.system('{0} --index {1}'.format(fetch_call, self.extracted[filename]))
        return True

    @Manager.milestone
    def download_database_files(self):
        filename = 'database_files/'
        result = self.download(filename)
        if result:
            destiny = self.downloaded[filename]
            path_pattern = '{0}{1[url]}{1[path]}Pfam{2}/{3}'
            server_tree = path_pattern.format(select_proto, self.manager.server,
                                         self.version, filename)
            server_tree = server_tree.split('//')[-1]
            origin = destiny + server_tree
            os.system('mv {:} {:}'.format(origin[:-1], destiny[:-1]))
            old_tree = destiny + server_tree.split('/')[0]
            shutil.rmtree(old_tree)
        return result


    @Manager.milestone
    def load_database(self):
        import loader
        loader.init_db(self)
        return True

    def upgrade(self):
        config = self.manager.config
        config["actual_version"] = self.version
        self.manager.config = config

    def prepare(self):
        exported = False
        while not exported:
            self.download_pfama_full()
            exported = self.export()
        self.clean_accession_numbers()
        self.reindex()
        self.download_database_files()
        self.load_database()

    def remove(self):
        # TODO: Should remove the database and the files of that specific
        # version.
        pass


from datetime import datetime, timedelta


# @app.run_every("day", "21:07")
dt = datetime.now() + timedelta(minutes=1)
@app.run_every("day", dt.strftime("%H:%M"))
def check_for_new_version():
    Manager().update()
