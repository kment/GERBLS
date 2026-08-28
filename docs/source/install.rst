Installing GERBLS
=================

.. attention::
    GERBLS currently requires **Python version 3.9** or above as well as a C++ compiler (explained
    below).

In addition to Python code, GERBLS also includes C++ code that is exposed to the Python interface
through `Cython`_. Regardless of the installation source below, the code is distributed as a source
distribution. **This means a C++ compiler (e.g., gcc, clang) is required to install GERBLS.** This
should not be an issue with most operating systems, but it bears mentioning nonetheless. GERBLS is
packaged using `setuptools`_, and building the package minimally requires `Cython`_ and `NumPy`_. 
These build-time dependencies will be set up automatically during the installation process.

.. _Cython: https://cython.org/
.. _NumPy: https://numpy.org/
.. _setuptools: https://setuptools.pypa.io/en/latest/

Installing from PyPI
--------------------
The easiest way to install GERBLS is by using `pip`. GERBLS comes in two versions, depending on your
needs.

Option 1: Basic install (sufficient for BLS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This version is sufficient to generate the BLS spectra, but some of the additional features (such as
fitting for a limb-darkened transit model) will not work without installing additional dependencies.
The only required dependencies are `NumPy`_ and `SciPy`_ (as well as a build-time dependency
on `Cython`_, as mentioned above). The basic version is installed by running: ::

    pip install gerbls

.. versionchanged:: 0.8
    SciPy is now a requirement for a minimal install.

.. _SciPy: https://scipy.org/

.. _full-install-label:

Option 2: Full install (with extras)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. warning::
    Future updates to dependent packages may cause GERBLS to break. If you encounter an error, try 
    reverting to a basic install and/or report the issue on the `GitHub repository`_.

The full install includes a small number of additional dependencies (currently `batman`_ and
`pytest`_) that enable some extra convenience features of GERBLS. To do a full install, simply run: ::

    pip install "gerbls[extras]"

.. _GitHub repository: https://github.com/kment/GERBLS/issues
.. _batman: http://lkreidberg.github.io/batman
.. _pytest: https://docs.pytest.org/en/stable/

Installing from source
----------------------
Alternatively, the GERBLS source repository can be cloned from GitHub and installed directly: ::

    git clone https://github.com/kment/GERBLS.git
    cd GERBLS
    pip install ".[extras]"

The ``[extras]`` option can be omitted for a basic install, as described above. A source install
also includes unit tests that can be run via `pytest`_ (which is included with ``[extras]``): ::

    pytest

Common issues
-------------

**1. I get an error stating "Python.h: No such file or directory"**

This error likely arises from the Python development headers not being installed on your system. The
solution is to install the ``python-dev`` package. A short tutorial can be found `here`_.

.. _here: https://betterstack.com/community/questions/how-to-fix-python-h-no-such-file-or-directory/