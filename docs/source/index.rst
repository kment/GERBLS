GERBLS documentation
====================

.. important::
   This documentation is currently undergoing improvements. Some package features be poorly
   documented, which will be improved in the near future.

**GERBLS** (Greatly Expedited Robust Box Least Squares) is a lightweight
fast-folding implementation of the BLS (Box Least Squares) algorithm. It is designed to facilitate
transiting planet searches in photometric data via an easy setup and fast runtimes.

GERBLS can outperform popular brute-force BLS implementations such as
`astropy.timeseries.BoxLeastSquares` by **over 10-20x** in runtime speed.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   getting_started
   basic_usage
   full_usage
   examples/fastbls
   api/index

* :ref:`genindex`