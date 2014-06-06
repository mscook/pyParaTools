pyParaTools - For exploring paramagnetic observables using theoretical models
=============================================================================

**I am no longer actively working on this project.**

Copyright (C) 2010-2014 Mitchell J Stanton-Cook

**Author:** Mitchell J Stanton-Cook
**Start Date:** 01/10
**Links:** http://comp-bio.anu.edu.au/mscook/PPT/
**Contact:** m.stantoncook@gmail.com

Please note you can cite the pyParaTools code using:

.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.10313.png
   :target: http://dx.doi.org/10.5281/zenodo.10313
   :alt: DOI

An alternate citation is::

    M. Stanton-Cook, X.-C. Su, G. Otting, T. Huber, http://compbio.anu.edu.au/mscook/PPT/


**Code from pyParaTools has been in incorporated into Xplor-NIH_ as of
version 2.32.**


Summary
-------

pyParaTools is a python module developed to work with paramagnetic nuclear
magnetic (NMR) observables in a more friendly manner.

The current version supports Pseudocontact shifts (PCS) Paramagnetic
Relaxation enhancement and Residual Dipolar Couplings.

pyParaTools provides a datastore for such assigned experimental data. It
can be used to calculate the expected expermintally determined values
from known/assumed parameters. It can be used to determine parameters
using non-linear least square fitting. Additional functions include
data exploration and utilities.

pyParaTools is easily extended. We ask that all modifications to pyParaTools
respect the licencing conditions. We would also like to hear how pyParaTools
has been used/modified.


Requirements
------------

You'll need to install (via pip, yum apt-get):
    * bio.pdb_,
    * Numpy_, and
    * Scipy_


References
----------

Please see::

    1) Bertini I, Luchinat C, Parigi G (2002) Magnetic susceptibility in
    paramagnetic NMR. Prog NMR Spectrosc 40:249–273

    2) Schmitz C, Stanton-Cook MJ, Su XC, Otting G, Huber T (2008)
    Numbat: an interactive software tool for fitting Δχ-tensors
    to molecular coordinates using pseudocontact shifts.
    J Biomol NMR 41:179–189

    3) Valafar H, Prestegard JH (2004) REDCAT: a residual dipolar coupling analysis
    tool. J Magn Reson 167:228–241

    4) Banci L, Bertini I, Cavallaro G, Giachetti A, Luchinat C, Parigi G (2004)
    Paramagnetism-based restraints for Xplor-NIH. J Biomol NMR 28:249–261


Licencing
---------

pyParaTools is liscenced under the `Educational Community License, Version 2.0`_
(ECL-2.0)

pyParaTools - For exploring paramagnetic observables using theoretical models 
pyParaTools  Copyright (C) 2010-2014  Mitchell J Stanton-Cook

pyParaTools comes with ABSOLUTELY NO WARRANTY; for details read LICENCE.txt
This is free software, and you are welcome to redistribute it
under certain conditions; read LICENCE.txt for details.



.. _Educational Community License, Version 2.0: http://opensource.org/licenses/ECL-2.0
.. _bio.pdb: http://www.biopython.org or from apt-get or yum
.. _Numpy: http://numpy.scipy.org or from apt-get or yum
.. _Scipy: http://www.scipy.org or from apt-get or yum
