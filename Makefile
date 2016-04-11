PYPREFIX_PATH=/usr
PYTHONLIBS=LD_LIBRARY_PATH=/usr/lib PATH=/usr/local/mysql/bin/:$(PATH)
PYTHONPATH=$(PYPREFIX_PATH)/bin/python
FIRST_EASYINSTALL=$(PYTHONLIBS) easy_install
PIP=pip
PYTHON=bin/python
EASYINSTALL=bin/easy_install
VIRTUALENV=virtualenv
SOURCE_ACTIVATE=$(PYTHONLIBS) . bin/activate;

unattended:
	@ (sudo ls 2>&1) >> tracking.log

ubuntu:
	@ sudo apt-get -y install zlibc libssl1.0.0 libbz2-dev libxslt1-dev libxml2-dev python-gevent libcurl4-openssl-dev wget mysql-server
	@ echo "[ assume       ] ubuntu distribution"

osx:
	@ sudo port install libevent mysql56-server
	@ echo "[ assume       ] osx system"

virtualenv:
	@ echo "[ installing   ] $(VIRTUALENV)"
	@ sudo $(FIRST_EASYINSTALL) virtualenv

bin/activate: requirements.txt
	@ echo "[ using        ] $(PYTHONPATH)"
	@ echo "[ creating     ] $(VIRTUALENV) with no site packages"
	@ ($(PYTHONLIBS) $(VIRTUALENV) --python=$(PYTHONPATH) --no-site-packages . 2>&1) >> tracking.log
	@ echo "[ installing   ] $(PIP) inside $(VIRTUALENV)"
	@ ($(SOURCE_ACTIVATE) $(EASYINSTALL) pip 2>&1) >> tracking.log
	@ echo "[ installing   ] $(PIP) requirements"
	@ $(SOURCE_ACTIVATE) $(PIP) install --upgrade pip
	@ $(SOURCE_ACTIVATE) $(PIP) install --upgrade distribute
	@ $(SOURCE_ACTIVATE) $(PIP) install -e  .
	@ $(SOURCE_ACTIVATE) $(PIP) install --default-timeout=100 -r requirements.development.txt
	@ touch bin/activate


deploy: bin/activate
	@ echo "[ deployed     ] the system was completly deployed"

show-version:
	@ $(SOURCE_ACTIVATE) $(PYTHON) --version

run:
	$(SOURCE_ACTIVATE) $(PYTHON) pfamserver

test:
	$(SOURCE_ACTIVATE) $(PYTHON) tests
	@ echo "[ tested       ] the system was completly tested"

shell:
	@ $(SOURCE_ACTIVATE) ipython
	@ echo "[ shell        ] the system was loaded into an ipython shell"

test-coverage-travis-ci:
	@ $(SOURCE_ACTIVATE) coverage run --source='pfamserver/' tests/__main__.py

test-coveralls:
	@ $(SOURCE_ACTIVATE) coveralls

test-coverage: test-coverage-travis-ci test-coveralls

pypi-register: test
	@ echo "[ record       ] package to pypi servers"
	@ ($(SOURCE_ACTIVATE) $(PYTHON) setup.py register -r pypi 2>&1) >> tracking.log
	@ echo "[ registered   ] the new version was successfully registered"

pypi-upload: test
	@ echo "[ uploading    ] package to pypi servers"
	@ ($(SOURCE_ACTIVATE) $(PYTHON) setup.py sdist upload -r https://pypi.python.org/pypi 2>&1) >> tracking.log
	@ echo "[ uploaded     ] the new version was successfully uploaded"

pypitest-register: test
	@ echo "[ record       ] package to pypi servers"
	@ $(SOURCE_ACTIVATE) $(PYTHON) setup.py register -r testpypi
	@ echo "[ registered   ] the new version was successfully registered"

pypitest-upload: test
	@ echo "[ uploading    ] package to pypi servers"
	$(SOURCE_ACTIVATE) $(PYTHON) setup.py sdist upload -r https://testpypi.python.org/pypi
	@ echo "[ uploaded     ] the new version was successfully uploaded"

clean:
	@ echo "[ cleaning     ] remove deployment generated files that doesn't exists in the git repository"
	@ rm -rf MANIFEST virtualenv* bin/ lib/ lib64 include/ build/ share setuptools-*.tar.gz get-pip.py tracking.log subversion .Python

create_boot_script:
	echo '#!/bin/sh\ncd $(shell pwd -P) ; $(SOURCE_ACTIVATE) $(PYTHON) -m pfamserver &' > /etc/init.d/S99pfamserver
	chmod +x /etc/init.d/S99pfamserver
	ln -s /etc/init.d/S99pfamserver /etc/rcS.d/
