.. _user_manual:

###########
User Manual
###########

***********
Quick Start
***********

This stub repo contains hooks for writing Abaqus :cite:`ABAQUS2019` subroutines, like those found in the `Abaqus UMAT
documentation`_, and a template UMAT c++ interface. However, this template repository does not yet have a meaningful c++
constitutive model to be the subject of a user manual.

This project is built and deployed to the `W-13 Python Environments`_ with continuous integration (CI) and continuous
deployment (CD). Most users will not need to build and install this project from source. Outside of the `W-13 Python
Environments`_, users may need to build and install directly from source. In that case, users are directed to the
:ref:`build` instructions.

With the `W-13 Python Environments`_, this project is installed in the Conda environment ``lib64`` and ``include``
directories, e.g. ``/path/to/my/conda/environment/{lib64,include}``. The template UMAT can be used with the following
Abaqus options after setting the system environment variable ``LD_LIBRARY_PATH``.

.. code:: bash

   $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path/to/conda/environment/lib64
   $ abaqus -job <my_input_file> -user path/to/conda/environment/lib64/cpp_stub_umat.o

Where the appropriate path can be found with

.. code:: bash

   $ find path/to/conda/environment -name "libcpp_stub.so"

For instance, with the W-13 "release" environment on ``sstelmo``

.. code:: bash

   $ find /projects/python/release -name "libcpp_stub.so"
   /projects/python/release/lib64/libcpp_stub.so
   $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/projects/python/release/lib64
   $ abaqus -job <my_input_file> -user /projects/python/release/lib64/cpp_stub_umat.o

As a convenience, the following code may be used to determine the correct, active Conda environment at Abaqus execution.
The following bash code is provided as an example for end users and not supported by this project. End users who wish to
learn more about bash scripting are directed to the online Bash documentation.

.. code:: bash

   # Get current conda environment information
   conda_env_path=$(conda info | grep "active env location" | cut -f 2 -d :)
   # Export the conda environment library path
   $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${conda_env_path}/lib64
   # Execute Abaqus with current Conda environment's installation of this project
   $ abaqus -job <my_input_file> -user ${conda_env_path}/lib64/cpp_stub_umat.o

***************************
Use after build from source
***************************

The template UMAT can be used after build with the following Abaqus options

.. code:: bash

   $ abaqus -job <my_input_file> -user relative/path/to/cpp_stub/build/src/cpp/cpp_stub_umat.o

It is strongly recommended that anyone building from source make use of the CMake ``--install`` options in a local Conda
environment. It is also possible to install to more traditional system paths, but this may require significantly more
background reading in relevant system administration.

Unless the template repository and all upstream c++ libraries are built and installed to a common system path it is
recommended that the subroutines are left in the project build directory. However, it is possible to copy the shared
library files to any other directory provided the upstream projects ``{error,vector,stress,solver,constitutive}_tools``
are present in the build directory, e.g.
``cpp_stub/build/_deps/{error,vector,stress,solver,constitutive}_tools-build/``.

.. code:: bash

   $ pwd
   /path/to/my/abaqus/job
   $ cp /path/to/cpp_stub/build/src/cpp/{cpp_stub_umat.o,libcpp_stub.so} .
   $ abaqus -job <my_input_file> -user cpp_stub_umat.o

******************************
Input File Material Definition
******************************

.. warning::

   Constitutive modeler health warning! The integration tests use a ``STATEV`` and ``PROPS`` length of one as the
   "incorrect" lengths to check the thrown exceptions. If your real constitutive model actually using a length of one
   for either vector, the integration test expectation must be updated.

cpp_stub requires 2 material constants and 2 state variables. The c++ cpp_stub interface, material constants, and state
variables are described in the :ref:`sphinx_api`. The fixed expectations for the abaqus interface are defined in the
"Variables" section of the :ref:`sphinx_api` for :ref:`cpp_stub_source`. A complete discussion about the constants and their
meaning is not included here. Instead users are directed to calibrated material parameters found in cpp_stub entries in the
`Granta/MIMS`_ `Material Database`_ :cite:`MIMS`. Material parameter calibration sets should be availble for download
with the correct Abaqus input file formatting from MIMS.

The cpp_stub project contains abaqus integration tests for the cpp_stub abaqus interface. These tests perform actual abaqus
simulations using the same dummy parameters used for unit and integration testing of the cpp_stub c++ code. The cpp_stub
Abaqus input files used for integration testing can be found in the cpp_stub source code repository with the following bash
command

.. code:: bash

   $ pwd
   /path/to/my/cpp_stub
   $ find . -path ./build -prune -false -o -name "*.inp"
   ./src/abaqus/single_element_c3d8.inp

The material definition from an integration test input file is included below for reference

.. warning::

   The material constants used in this example material definition are *NOT* calibrated for any real material data.
   For calibrated material parameters, see the cpp_stub entry for materials found in the `Granta/MIMS`_ `Material
   Database`_ :cite:`MIMS`.

.. literalinclude:: ../../src/abaqus/single_element_c3d8.inp
   :linenos:
   :lines: 42-50
