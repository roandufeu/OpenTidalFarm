Installation
============

It is recommended to use OpenTidalFarm on the `Ubuntu`_ operating system.

Quick install
-------------

Manual installation (recommended)
*********************************

For the manual installation of OpenTidalFarm you need to have the Ubuntu Linux.
Then you need to install these OpenTidalFarm dependencies:

- `FEniCS`_ (Follow the Ubuntu PPA installation)
- `dolfin-adjoint`_ (Follow the Ubuntu PPA installation)
- `SciPy`_

   ``pip install scipy``

- `pyyaml`_

   ``pip install pyyaml``

- `Uptide`_

   ``pip install git+git://github.com/stephankramer/uptide.git``

- `UTM`_

   ``pip install utm``


Finally, you can install OpenTidalFarm with this command:

.. code-block:: bash

   pip install git+git://github.com/OpenTidalFarm/OpenTidalFarm.git

Once finished, you can test the installation with:

.. code-block:: bash

    python -c "import opentidalfarm"

If no errors occur, your installation was succesfull.

One can download the examples with all mesh data with these commands:

.. code-block:: bash

   git clone git@github.com:OpenTidalFarm/OpenTidalFarm.git
   cd OpenTidalFarm
   git submodule init
   git submodule update

The examples are then stored in the "OpenTidalFarm/examples" directory.

.. _Ubuntu: http://www.ubuntu.com/
.. _FEniCS: http://fenicsproject.org/download/
.. _dolfin-adjoint: http://www.dolfin-adjoint.org/download/index.html
.. _SciPy >=0.11: https://github.com/scipy/scipy
.. _Uptide: https://github.com/stephankramer/uptide
.. _UTM: https://pypi.python.org/pypi/utm
.. _Download OpenTidalFarm: https://github.com/funsim/OpenTidalFarm/zipball/master
.. _Issue tracker: https://github.com/OpenTidalFarm/OpenTidalFarm/issues
.. _SciPy: http://www.scipy.org
.. _pyyaml: http://pyyaml.org

Automatic installation using hashdist (experimental)
****************************************************

Following command will install OpenTidalFarm and all its dependencies into a isolated environment on your computer.

.. code-block:: bash

   curl -s https://bitbucket.org/simon_funke/fenics-developer-tools/raw/master/opentidalfarm-install.sh | bash

Select the option "[1] latest stable version of FEniCS" during the installation.

Once finished, you can test the installation with:

.. code-block:: bash

    python -c "import opentidalfarm"

If no errors occur, your installation was succesfull.

If you have any problems with the installation, please use our `Issue tracker`_.
