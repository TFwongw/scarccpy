========================
Installation of SCARCCpy
========================

Use of Python3.10 is required.

Use of `Python virtual environment <https://realpython.com/python-virtual-environments-a-primer>`_ is encouraged.

To set up in Minnesota Supercomputing Institute(MSI), run the following commands:


.. code-block:: console
    
    $ module load python
    $ python3 -m ensurepip
    $ python3 -m pip install virtualenv
    $ python3 -m venv scarccpy

.. note::

   Do not use `module load python3` in MSI. The conda environment only supports up to Python 3.8. Whereas SCARCCpy requires Python 3.10.

For usage instructions, please see the `documentation <https://scarccpy.readthedocs.io/en/latest/>`_.