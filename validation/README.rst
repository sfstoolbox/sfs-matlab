Validation scripts
==================

This folder contains test scripts for the most important functions of the SFS
Toolbox.

In order to check if no error during execution occurs, run:

.. sourcecode:: Matlab

    >> test_all(0)

which should return a ``1`` if everything runs smoothly.

In order to get a detailed overview if the algorithms return desired results you
have to check the resulting plots via inspection. To do so run

.. sourcecode:: Matlab

    >> test_all(1)

and follow the instructions presented in the command prompt.
