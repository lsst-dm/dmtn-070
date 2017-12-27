.. image:: https://img.shields.io/badge/dmtn--069-lsst.io-brightgreen.svg
   :target: https://dmtn-069.lsst.io
.. image:: https://travis-ci.org/lsst-dm/dmtn-069.svg
   :target: https://travis-ci.org/lsst-dm/dmtn-069
..
  Uncomment this section and modify the DOI strings to include a Zenodo DOI badge in the README
  .. image:: https://zenodo.org/badge/doi/10.5281/zenodo.#####.svg
     :target: http://dx.doi.org/10.5281/zenodo.#####

#################################################
Report on Summer 2014 Production: Analysis of DCR
#################################################

DMTN-069
========

The goals of this Summer 2014 (S14) task were to understand the scope
of the differential chromatic refraction (DCR) issue using a realistic
range of stellar spectral energy distributions (SEDs).  We used LSST
catSim's all--sky catalog of stellar sources, including their number
counts and magnitude distributions, to estimate how many sources will
be detected at 5--sigma in the LSST $ugri$--bands.  The SED of each of
these sources was used to model the per--passband refraction, and then
the differential refraction with respect to a single reference SED, as
a function of airmass.

**Links:**

- Publication URL: https://dmtn-069.lsst.io
- Alternative editions: https://dmtn-069.lsst.io/v
- GitHub repository: https://github.com/lsst-dm/dmtn-069
- Build system: https://travis-ci.org/lsst-dm/dmtn-069


Build this technical note
=========================

You can clone this repository and build the technote locally with Latex.
You must have `lsst-texmf`_ installed.

.. code-block:: bash

   git clone https://github.com/lsst-dm/dmtn-069
   cd dmtn-069
   make


Editing this technical note
===========================

You can edit the ``DMTN-069.tex`` file, which is a latex document.

Remember that images and other types of assets should be stored in the ``_static/`` directory of this repository.
See ``_static/README.rst`` for more information.

The published technote at https://dmtn-069.lsst.io will be automatically rebuilt whenever you push your changes to the ``master`` branch on `GitHub <https://github.com/lsst-dm/dmtn-069>`_.

****

Copyright 2017 University of Washington

This work is licensed under the Creative Commons Attribution 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.

.. _this repo: ./DMTN-069.tex
.. _lsst-texmf: https://lsst-texmf.lsst.io
