config = {}
config['DEBUG'] = True
config['SECRET_KEY'] = 'some_secret_key!'
config.update(
    SQLALCHEMY_DATABASE_URI = 'mysql+pymysql://root:root@192.168.101.100:3306/pfamserver',
    #SQLALCHEMY_DATABASE_URI = 'mysql+pymysql://root:root@localhost:3306/pfamserver',
    #SQLALCHEMY_DATABASE_URI = 'sqlite:///test.db',
    #SQLALCHEMY_DATABASE_URI = 'postgresql://postgres:postgres@localhost:5432/pfamserver',
    SESSION_TYPE = 'sqlalchemy',
)
