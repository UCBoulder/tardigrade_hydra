.. targets-start-do-not-remove

.. _`AEA Conda channel`: https://aea.re-pages.lanl.gov/developer-operations/aea_compute_environment/aea_compute_environment.html#aea-conda-channel
.. _`AEA compute environment`: https://aea.re-pages.lanl.gov/developer-operations/aea_compute_environment/aea_compute_environment.html#
.. _Anaconda Documentation: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html
.. _BOOST: https://www.boost.org/doc/libs/1_53_0/
.. _`Conda`: https://docs.conda.io/en/latest/
.. _CMake: https://cmake.org/cmake/help/v3.14/
.. _CMake add_custom_target: https://cmake.org/cmake/help/latest/command/add_custom_target.html
.. _CMake fetch_content: https://cmake.org/cmake/help/latest/module/FetchContent.html
.. _Doxygen: https://www.doxygen.nl/manual/docblocks.html
.. _Eigen: https://eigen.tuxfamily.org/dox/
.. _Sphinx: https://www.sphinx-doc.org/en/master/
.. _Breathe: https://breathe.readthedocs.io/en/latest/
.. _PEP-8: https://www.python.org/dev/peps/pep-0008/
.. _pipreqs: https://github.com/bndr/pipreqs
.. _LaTeX: https://www.latex-project.org/help/documentation/
.. _upstream repository: https://re-git.lanl.gov/aea/stub-repositories/tardigrade-hydra
.. _Material Models: https://re-git.lanl.gov/aea/material-models
.. _UNIX group: https://ddw-confluence.lanl.gov/pages/viewpage.action?pageId=150929410

.. targets-end-do-not-remove

################
Tardigrade Hydra
################

*******************
Project Description
*******************

.. project-brief-start-do-not-remove

A C++ framework to develop finite deformation material models

.. project-brief-end-do-not-remove

Information
===========

* Documentation (``main`` branch): https://aea.re-pages.lanl.gov/stub-repositories/tardigrade-hydra/

* Wiki: https://re-git.lanl.gov/aea/stub-repositories/tardigrade-hydra/-/wikis/home

Developers
==========

* Nathan Miller: nathanm@lanl.gov

************
Gitlab CI/CD
************

    **NOTE**

    The repository setup has moved out of the README and into the HTML
    documentation. You can find the Gitlab project setup guide here:
    https://aea.re-pages.lanl.gov/stub-repositories/tardigrade-hydra/gitlab_setup.html

************
Dependencies
************

.. dependencies-start-do-not-remove

For convenience, the minimal Conda environment requirements for project development are included in ``environment.txt``.
A minimal anaconda environment for building the documentation can be created from an existing anaconda installation with
the following commands.

.. code-block:: bash

   $ conda create --name tardigrade-hydra-env --file environment.txt --channel file:///projects/aea_compute/aea-conda

You can learn more about Anaconda Python environment creation and management in
the `Anaconda Documentation`_.

The build, test, and run time requirements are subsets of the development environment requirements found in
``environment.txt``. This project builds and deploys as a Conda package to the `AEA Conda channel`_. The Conda recipe
and build, test, run time requirements are found in the ``recipe/`` directory.

.. dependencies-end-do-not-remove

**************
Build and Test
**************

.. build-start-do-not-remove

This project is built with `CMake`_ and uses `Sphinx`_ to build the
documentation with `Doxygen`_ + `Breathe`_ for the c++ API.

Environment variables
=====================

This project's `CMake`_ configuration accepts two build type strings: 'Release' and 'conda-test'. The first is used
during the Gitlab-CI ``fast-test`` job to ensure that the project uses installed libraries correctly. The latter is used
during the Gitlab-CI ``conda-build`` job to limit the test phase to the as-installed project files.

The build type can be set with the ``-DCMAKE_BUILD_TYPE=<build type string>`` during project configuration. Both build
types will require the upstream dependent libraries

* ``abaqus_tools``: https://re-git.lanl.gov/aea/material-models/abaqus_tools
* ``error_tools``: https://re-git.lanl.gov/aea/material-models/error_tools

to be installed and found in the user's environment. If the build type string doesn't match those previously listed, the
CMake project will build missing upstream libraries with the `CMake fetch_content`_ feature. The 'conda-test' build type
excludes the project libraries from the build configuration and will attempt to find the project libraries in the user's
environment to perform the project unit and integration tests against the as-installed project files.

Build on sstelmo
================

1) Activate the correct python environment

   .. code-block:: bash

      $ module use /projects/aea_compute/modulefiles
      $ module load tardigrade-hydra-env

2) Create a build directory

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade-hydra/

      $ mkdir build
      $ cd build

3) Configure ``cmake``

       This step only needs to be performed once unless you need to specify a new CMake configuration for a re-build.
       Most command line arguments and environment variables are stored in the CMake cache. Anything found in cache will
       not be re-configured unless you build the ``rebuild_cache`` target or clobber the build directory.

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade-hydra/build
      $ cmake ..

4) Display target options

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade-hydra/build
      $ cmake --build . --target help

4) Build various portions of the project

       Most of the project will re-build only as necessary after source updates. Some portions of the documentation
       require a ``cmake --build . --target clean`` after documentation source file updates to force a re-build.

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade-hydra/build

      # Build everything (either or)
      $ cmake --build .
      $ cmake --build . --target all

      # Build the c++ primary libraries by target name(s)
      $ cmake --build . --target tardigrade-hydra tardigrade-hydra_umat

      # Build the c++ primary libraries by sub-directory
      $ cmake --build src/cpp

5) Locate build files

       The build directory structure may change between version releases. Developers and users are encouraged to become
       familiar with the bash ``find``, ``grep``, and ``tree`` commands to locate build files.

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade-hydra/build

      # find c++ libraries and ignore intermediate files with similar extensions
      $ find . \( -name "*.o" -o -name "*.so" -o -name "*.a" \) | grep -vE "\.cpp\."

6) Clean build directory to force a re-build

       **HEALTH WARNING**

       The abaqus input files and bash scripts used for integration testing are
       built with the `CMake add_custom_target`_ feature. Consequently, the integration
       test target is *always considered out of date*. The integration test target
       copies all registered input files and the integration test bash script from
       source to build directory. This means the file copy operation is always
       performed when the integration test target is requested in the cmake build
       command, e.g. ``cmake --build .`` or ``cmake --build src/abaqus/tests``. This
       operation is computationally inexpensive with respect to building the
       ``tardigrade-hydra`` source code.

       Input files are registered in the ``src/abaqus/tests/CMakeLists.txt`` file
       under the ``ABAQUS_INPUT_FILES`` CMake variable.

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade-hydra/build

      $ cmake --build . --target clean

Test on sstelmo
===============

4) Build tests of the project

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade-hydra/build

      # Build c++ and abaqus tests by target name(s)
      $ cmake --build . --target test_tardigrade-hydra test_abaqus_integration

      # Build c++ and abaqus tests by sub-directories
      $ cmake --build src/cpp/tests
      $ cmake --build src/abaqus/tests

5) Run the tests

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade-hydra/build

      # Run ctest
      $ ctest

      # Results print to screen
      # View details of most recent test execution including failure messages
      $ less Testing/Temporary/LastTest.log

Building the documentation
==========================

    **HEALTH WARNING**

    The sphinx API docs are a work-in-progress. The doxygen API is much more useful.

    * Documentation (``main`` branch): https://aea.re-pages.lanl.gov/stub-repositories/tardigrade-hydra/doxygen

To build just the documentation pick up the steps here:

2) Create the build directory and move there

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade-hydra/
      $ mkdir build/
      $ cd build/

3) Run cmake configuration

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade-hydra/build/
      $ cmake ..

4) Build the docs

   .. code-block:: bash

      $ cmake --build . --target Sphinx

5) Documentation builds to:

   .. code-block:: bash

      tardigrade-hydra/build/docs/sphinx/html/index.html

6) Display docs

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade-hydra/build/
      $ firefox docs/sphinx/html/index.html &

7) While the Sphinx API is still a WIP, try the doxygen API

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade-hydra/build/
      $ firefox docs/doxygen/html/index.html &

*******************
Install the library
*******************

Build the entire before performing the installation.

4) Build the entire project

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade-hydra/build
      $ cmake --build .

5) Install the library

   .. code-block:: bash

      $ pwd
      /path/to/tardigrade-hydra/build
      $ cmake --install . --prefix path/to/root/install

      # Example local user (non-admin) Linux install
      $ cmake --install . --prefix /home/$USER/.local

      # Example install to conda environment
      $ conda activate my_env
      $ cmake --install . --prefix ${CONDA_PREFIX}

.. build-end-do-not-remove

***********************
Contribution Guidelines
***********************

.. contribution-start-do-not-remove

Git Commit Message
==================

Begin Git commit messages with one of the following headings:

* BUG: bug fix
* DOC: documentation
* FEAT: feature
* MAINT: maintenance
* TST: tests
* REL: release
* WIP: work-in-progress

For example:

.. code-block:: bash

   git commit -m "DOC: adds documentation for feature"

Git Branch Names
================

When creating branches use one of the following naming conventions. When in
doubt use ``feature/<description>``.

* ``bugfix/\<description>``
* ``feature/\<description>``
* ``release/\<description>``

reStructured Text
=================

`Sphinx`_ reads in docstrings and other special portions of the code as reStructured text. Developers should follow
styles in this `Sphinx style guide
<https://documentation-style-guide-sphinx.readthedocs.io/en/latest/style-guide.html#>`_.

Style Guide
===========

This project does not yet have a full style guide. Generally, wherever a style
can't be inferred from surrounding code this project falls back to `PEP-8`_-like
styles. There are two notable exceptions to the notional PEP-8 fall back:

1. `Doxygen`_ style docstrings are required for automated, API from source documentation.
2. This project prefers expansive whitespace surrounding parentheses, braces, and
   brackets.

   * No leading space between a function and the argument list.
   * One space following an open paranthesis ``(``, brace ``{``, or bracket
     ``[``
   * One space leading a close paranthesis ``)``, brace ``}``, or bracket ``]``

An example of the whitespace style:

.. code-block:: bash

   my_function( arg1, { arg2, arg3 }, arg4 );

The following ``sed`` commands may be useful for updating white space, but must
be used with care. The developer is recommended to use a unique git commit
between each command with a corresponding review of the changes and a unit test
run.

* Trailing space for open paren/brace/bracket

  .. code-block:: bash

     sed -i 's/\([({[]\)\([^ ]\)/\1 \2/g' <list of files to update>

* Leading space for close paren/brace/bracket

  .. code-block:: bash

     sed -i 's/\([^ ]\)\([)}\]]\)/\1 \2/g' <list of files to update>

* White space between adjacent paren/brace/bracket

  .. code-block:: bash

     sed -i 's/\([)}\]]\)\([)}\]]\)/\1 \2/g' <list of files to update>

.. contribution-end-do-not-remove
