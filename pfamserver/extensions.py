from __future__ import unicode_literals

from flask_compress import Compress
from future.standard_library import install_aliases
install_aliases()

from flask_caching import Cache
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
from flask_bootstrap import Bootstrap
from flask_mail import Mail
from flask_webpack import Webpack
from flask_wtf.csrf import CsrfProtect
import requests
from six.moves.urllib.parse import quote
from base64 import b64decode
import json


class PfamServerExtension:
    host = None
    headers = {'X-Requested-With': 'XMLHttpRequest'}

    def __init__(self):
        pass

    def init_app(self, app):
        self.host = app.config['PFAMSERVER_HOST']

    def obtain(self, partial):
        request = requests.get('{host}/api/{partial}'.format(host=self.host, partial=partial),
                                headers=self.headers)

        request.encoding = 'UTF-8'
        return request

    def pfams_from_uniprot(self, uniprot):
        partial = 'query/pfam_uniprot/{uniprot}'.format(uniprot=quote(uniprot))
        res = self.obtain(partial)
        if res.status_code == requests.codes.ok and res.text:
            data = res.json()
            pfams = []
            if 'output' in data:
                pfams = data["output"]
                for pfam in pfams:
                    pfam[u'id'] = u'{0[pfamA_acc]}/{0[seq_start]}-{0[seq_end]}'.format(pfam)
            return pfams

    def queryTable(self, table, column, value):
        filtered = json.dumps({'single': True,
                               'filters': [{'name': column, 'op': 'eq', 'val': value}]},
                              ensure_ascii=False).encode('utf8')
        return self.obtain('{table}?q={filtered}'.format(table=table, filtered=filtered))

    def uniprot(self, id=None, acc=None):
        if (id and acc) or not (id or acc):
            return None
        if acc:
            res = self.queryTable('Uniprot', 'uniprot_acc', acc)
        else:
            res = self.queryTable('Uniprot', 'uniprot_id', id)
        if res.status_code == requests.codes.ok and res.text:
            return res.json()

    def pfam(self, id=None, acc=None):
        if (id and acc) or not (id or acc):
            return ''
        if acc:
            res = self.queryTable('PfamA', 'pfamA_acc', acc)
        else:
            res = self.queryTable('PfamA', 'pfamA_id', id)
        if res.status_code == requests.codes.ok and res.text:
            return res.json()

    def stockholm_from_pfam(self, pfam):
        partial = 'query/stockholm_pfam/{pfam}'.format(pfam=quote(pfam))
        res = self.obtain(partial)
        if res.status_code == requests.codes.ok and res.text:
            data = res.json()
            if data['output']:
                return b64decode(data["output"])

    def sequencedescription_from_pfam(self, pfam, with_pdb=True):
        partial = 'query/sequencedescription_pfam/{pfam}?with_pdb={with_pdb}'.format(pfam=pfam,
                                                                                     with_pdb=with_pdb)
        res = self.obtain(partial)
        if res.status_code == requests.codes.ok:
            return res.json()["output"]

    def pdbs_from_sequencedescription(self, sequencedescription):
        sequencedescription = sequencedescription.replace('/', ',').replace('-', ',')
        partial = 'query/pdb_sequencedescription/{sequencedescription}'.format(sequencedescription=sequencedescription)
        res = self.obtain(partial)
        if res.status_code == requests.codes.ok:
            pdbs = res.json()["output"]
            for pdb in pdbs:
                pdb[u'id'] = u'{0[pdb_id]}/{0[chain]}/{0[pdb_res_start]}'.format(pdb)
            return pdbs


class PDBServerExtension:
    def __init__(self):
        pass

    def init_app(self, app):
        pass

    def file_from_pdb_id(self, pdb_id):
        url = 'http://files.rcsb.org/download/{pdb_id}.pdb'.format(pdb_id = pdb_id)
        return requests.get(url).content


bootstrap = Bootstrap()
cache = Cache(config={'CACHE_TYPE': 'simple'})
db = SQLAlchemy()
compress = Compress()
mail = Mail()
migrate = Migrate()
csrf = CsrfProtect()
webpack = Webpack()
pfamserver = PfamServerExtension()
pdbserver = PDBServerExtension()
