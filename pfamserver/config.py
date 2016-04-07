from distutils.sysconfig import get_python_lib
from os import getenv


config = {}
config['DEBUG'] = eval(getenv('DEBUG', 'True'))
config['SECRET_KEY'] = 'some_secret_key!'
config['LIB_PATH'] = get_python_lib()
config['ROOT_PATH'] = "{:s}/pfam_data/".format(config['LIB_PATH'])
config.update(
    #SQLALCHEMY_DATABASE_URI = 'mysql+pymysql://root:root@192.168.101.100:3306/pfamserver',
    SQLALCHEMY_DATABASE_URI = 'mysql+pymysql://{:}'.format(
        getenv('DB', 'root:root@localhost:3306/pfamserver')),
    #SQLALCHEMY_DATABASE_URI = 'sqlite:///test.db',
    #SQLALCHEMY_DATABASE_URI = 'postgresql://postgres:postgres@localhost:5432/pfamserver',
    SESSION_TYPE = 'sqlalchemy',
)
