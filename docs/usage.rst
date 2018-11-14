=====
Usage
=====

To use MPC Utilities in a project::

    import mpcutilities

=======
Kepcart
=======

Details of the functions in the :mod:`kepcart` module for transforming between Keplerian and Cartesian coordinates.
Note that one of the functions, ``kepState2cartPV`` is currently commented-out and unusable.

.. automodule:: mpcutilities.kepcart
   :members:		

==================
Kepcart Structures
==================

Details of the c-types Structures defined in the :mod:`classes` module that are used for passing information back-and-forth between python & C:

.. automodule:: mpcutilities.classes
   :members:		

 
=========
PhysConst
=========

Details of the rotation functions which are currently located in :mod:`phys_const` module  (for some reason):

.. automodule:: mpcutilities.phys_const
   :members:		

=========
obs80
=========

Details of the functions for parsing 80-character-observations (written by Sonia Keys) and located in :mod:`obs80` module.
This also contains "heliocentric-annotation'' code that used to be in "obs80hc".

.. automodule:: mpcutilities.obs80
   :members:		

=========
ele220
=========

Details of the functions for parsing 220-character-orbits (written by Sonia Keys) and located in :mod:`ele220` module.

.. automodule:: mpcutilities.ele220
   :members:		



==============
leapsec module
==============

Archaic / Unused leap-second functionality: documentation included in case it needs to be tracked-down

.. automodule:: mpcutilities.leapsec


==============
obscode module
==============

Archaic / Unused obscode functionality: documentation included in case it needs to be tracked-down

.. automodule:: mpcutilities.obscode
