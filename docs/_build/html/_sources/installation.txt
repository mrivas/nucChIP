.. _installation:

******************************
Prerequisites and installation
******************************

nucChIP source code and installation scripts are available at

	`nucChIP github page`_

.. _nucChIP github page: https://github.com/mrivas/nucChIP


To use nucChIP, you need at least version 2.5 of Python_ (Python 3 does not work yet), 

.. _Python: http://www.python.org/

nucChIP installation package will take care of all dependencies. For users unfamiliar downloading source code from github or using python installation scripts, detailed instructions follow.


Installation on Linux
=====================

Download the *source* package compressed as zip file from `nucChIP github page`_, unzip the file, and go into the directory with the unzipped files::
	
	wget https://github.com/mrivas/nucChIP/archive/master.zip nucChIP.zip
	unzip nucChIP.zip
	cd nucChIP

To install nucChIP for the user currently logged in type::

	python setup.py install --user

To make nucChIP available to all users, use instead::

	python setup.py build
	sudo python setup.py install

To test the installation, change to another director than the build directory, start Python (by typing ``python``) and then try whether typing ``import nucChIP`` causes an error meesage.

