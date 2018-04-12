from distutils.sysconfig import get_python_lib
from os import getenv


config = {}
config['DEBUG'] = eval(getenv('DEBUG', 'True'))
config['SECRET_KEY'] = 'some_secret_key!'
config['LIB_PATH'] = get_python_lib()

# path to store and look for pfam data (+100GB)
config['ROOT_PATH'] = '/home/mistic/pfam_data/' # it is important to put the slash at the end..
config['DBMANUAL_PATH'] = '/home/mistic/pfam31/database_files'
config['TEMP'] = '/tmp/'
config['PFAMSCAN_PATH'] = '/home/mistic/pfam_data/31.0/pfamscan'
config.update(
    #SQLALCHEMY_DATABASE_URI = 'mysql+pymysql://root:root@192.168.101.100:3306/pfamserver',
    SQLALCHEMY_DATABASE_URI = 'mysql+pymysql://{:}'.format(
        getenv('DB', 'root:mistic2016@localhost:3306/pfamserver')),
    #SQLALCHEMY_DATABASE_URI = 'sqlite:///test.db',
    #SQLALCHEMY_DATABASE_URI = 'postgresql://postgres:postgres@localhost:5432/pfamserver',
    SESSION_TYPE = 'sqlalchemy',
)
