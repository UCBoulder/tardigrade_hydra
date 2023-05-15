.. _user_manual:

###########
User Manual
###########

***********
Quick Start
***********

This stub repo contains hooks for writing Abaqus :cite:`ABAQUS2022` subroutines, like those found in the `Abaqus UMAT
documentation`_, and a template UMAT c++ interface. However, this template repository does not yet have a meaningful c++
constitutive model to be the subject of a user manual.

This project is built and deployed to the `AEA Conda channel`_ with continuous integration (CI) and continuous
deployment (CD). The `AEA compute environment`_ installs this project from the `AEA Conda channel`_. Most users will not
need to build or install this project from source. Outside of the `AEA compute environment`_, users can install directly
from the `AEA Conda channel`_. In rare cases, users may need to build from source and are directed to the :ref:`build`
instructions.

From the `AEA Conda channel`_, this project is installed in the Conda environment ``lib64`` and ``include`` directories,
e.g. ``${CONDA_PREFIX}/{lib64,include}``. When the `AEA compute environment`_ module files are used for environment
activation, the template UMAT can be used with the following Abaqus options.

.. code:: bash

   $ abaqus -job <my_input_file> -user ${CONDA_PREFIX}/lib64/tardigrade-hydra_umat.o

Where the appropriate path can be confirmed with

.. code:: bash

   $ find ${CONDA_PREFIX} -name "libtardigrade-hydra.so"

For instance, with the "aea-release" environment on ``sstelmo``

.. code:: bash

   # See the AEA compute environment documentation to confirm the preferred activation command
   $ module use /projects/aea_compute/modulefiles
   $ module load aea-beta

   $ echo ${CONDA_PREFIX}
   /projects/aea_compute/aea-release
   $ find ${CONDA_PREFIX} -name "tardigrade-hydra_umat.o"
   /projects/aea_compute/aea-release/lib64/tardigrade-hydra_umat.o
   $ abaqus -job <my_input_file> -user ${CONDA_PREFIX}/lib64/tardigrade-hydra_umat.o

If the `AEA compute environment`_ module files are not used, the user must set their ``LD_LIBRARY_PATH`` manually. As a
convenience, the following code may be used to determine the active Conda environment at Abaqus execution. The following
bash code is provided as an example for end users and not supported by this project. End users who wish to learn more
about bash scripting are directed to the online Bash documentation.

.. code:: bash

   # Export the conda environment library path
   $ export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:${CONDA_PREFIX}/lib64:$LD_LIBRARY_PATH

   # Execute Abaqus with current Conda environment's installation of this project
   $ abaqus -job <my_input_file> -user ${CONDA_PREFIX}/lib64/tardigrade-hydra_umat.o

***************************
Use after build from source
***************************

The template UMAT can be used after build with the following Abaqus options

.. code:: bash

   $ abaqus -job <my_input_file> -user relative/path/to/tardigrade-hydra/build/src/cpp/tardigrade-hydra_umat.o

It is strongly recommended that anyone building from source make use of the CMake ``--install`` options in a local Conda
environment as in :ref:`build`. It is also possible to install to more traditional system paths, but this may require
significantly more background reading in relevant system administration.

Unless the template repository and all upstream c++ libraries are built and installed to a common system path it is
recommended that the subroutines are left in the project build directory. However, it is possible to copy the shared
library files to any other directory provided the upstream projects ``{error,vector,stress,solver,constitutive}_tools``
are present in the build directory, e.g.
``tardigrade-hydra/build/_deps/{error,vector,stress,solver,constitutive}_tools-build/``.

.. code:: bash

   $ pwd
   /path/to/my/abaqus/job
   $ cp /path/to/tardigrade-hydra/build/src/cpp/{tardigrade-hydra_umat.o,libtardigrade-hydra.so} .
   $ abaqus -job <my_input_file> -user tardigrade-hydra_umat.o

******************************
Input File Material Definition
******************************

.. warning::

   Constitutive modeler health warning! The integration tests use a ``STATEV`` and ``PROPS`` length of one as the
   "incorrect" lengths to check the thrown exceptions. If your real constitutive model actually using a length of one
   for either vector, the integration test expectation must be updated.

tardigrade-hydra requires 2 material constants and 2 state variables. The c++ tardigrade-hydra interface, material constants, and state
variables are described in the :ref:`sphinx_api`. The fixed expectations for the abaqus interface are defined in the
"Variables" section of the :ref:`sphinx_api` for :ref:`tardigrade-hydra_source`. A complete discussion about the constants and their
meaning is not included here. Instead users are directed to calibrated material parameters found in tardigrade-hydra entries in the
`Granta/MIMS`_ `Material Database`_ :cite:`MIMS`. Material parameter calibration sets should be availble for download
with the correct Abaqus input file formatting from MIMS.

The tardigrade-hydra project contains abaqus integration tests for the tardigrade-hydra abaqus interface. These tests perform actual abaqus
simulations using the same dummy parameters used for unit and integration testing of the tardigrade-hydra c++ code. The tardigrade-hydra
Abaqus input files used for integration testing can be found in the tardigrade-hydra source code repository with the following bash
command

.. code:: bash

   $ pwd
   /path/to/my/tardigrade-hydra
   $ find . -path ./build -prune -false -o -name "*.inp"
   ./src/abaqus/single_element_c3d8.inp

The material definition from an integration test input file is included below for reference

.. warning::

   The material constants used in this example material definition are *NOT* calibrated for any real material data.
   For calibrated material parameters, see the tardigrade-hydra entry for materials found in the `Granta/MIMS`_ `Material
   Database`_ :cite:`MIMS`.

.. literalinclude:: ../../src/abaqus/single_element_c3d8.inp
   :linenos:
   :lines: 42-50
