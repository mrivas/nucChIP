.. _installation:

******************************
Prerequisites and installation
******************************

nucChIP source code, installation scripts, and source code of this website are available at

	`nucChIP github page`_

.. _nucChIP github page: https://github.com/mrivas/nucChIP

Prerequisites
=============

To use nucChIP, you need at least version 2.5 of `Python <http://www.python.org/>`_ (Python 3 doesn't work yet).


Installation on Linux
=====================

nucChIP installation package will take care of all dependencies. You just need to download the *source* files, unzip them, and go into the unzipped directory.

.. code-block:: bash

   $ wget https://github.com/mrivas/nucChIP/archive/master.zip nucChIP.zip
   $ unzip nucChIP.zip
   $ cd nucChIP

To install nucChIP locally ( only for the current user in case you don't have super-user permission) type

.. code-block:: bash

   $ python setup.py install --user

To make nucChIP available to all users, use instead

.. code-block:: bash

   $ python setup.py build
   $ sudo python setup.py install

To test the installation

.. code-block:: bash

   $ python
   >> import nucChIP

if no errors are reported, your are all set.
