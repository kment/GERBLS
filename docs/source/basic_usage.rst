Basic usage
===========
Basic usage is described on the `GERBLS GitHub page`_. However, **running GERBLS this way has
limited functionality** and is generally not recommended, other than for general testing purposes.
**To leverage the full functionality of GERBLS,** refer to the fastbls_ example.

.. _fastbls: examples/fastbls.ipynb
.. _GERBLS GitHub page: https://github.com/kment/GERBLS

Basic usage is implemented via a wrapper function :func:`gerbls.run_bls` that returns a dictionary
containing the generated BLS spectrum:

.. autofunction:: gerbls.run_bls
    :no-index: