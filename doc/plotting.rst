.. _sec-plotting:

Plotting
========

The Toolbox provides you with a function for plotting your simulated sound
fields (``plot_sound_field()``) and adding loudspeaker symbols to the figure
(``draw_loudspeakers()``). If you have gnuplot installed, you can use the
functions ``gp_save_matrix()`` and ``gp_save_loudspeakers()`` to save your data
in a way that it can be used with gnuplot. An example use case can be found `at
this plot of a plane wave`_ which includes the Matlab/Octave code to generate
the data and the gnuplot script for plotting it.

.. _at this plot of a plane wave: https://github.com/hagenw/phd-thesis/tree/master/02_theory_of_sound_field_synthesis/fig2_04

.. vim: filetype=rst spell:
