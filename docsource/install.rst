Installation instructions
=========================

Installing PISCES from git
--------------------------------

.. code:: shell

  $ pip install -e git+https://github.com/Novartis/pisces.git#egg=novartis-pisces

Installing PISCES from PyPI
---------------------------

.. code:: shell

  $ pip install novartis-pisces  


Testing your PISCES installation
--------------------------------

.. code:: shell

	# run the python module unit tests
	$ python -m tox

.. note::

    PISCES will automatically use the system installation of Salmon if it is available in your PATH. If not found, it will fall back to the bundled version included with PISCES.

