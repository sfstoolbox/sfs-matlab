Sound Field Synthesis Toolbox for Matlab
========================================

The SFS Toolbox for Matlab gives you the possibility to play around with sound
field synthesis methods like wave field synthesis (WFS) or near-field
compensated higher order Ambisonics (NFC-HOA).  There are functions to simulate
monochromatic sound fields for different secondary source (loudspeaker) setups,
time snapshots of full band impulses emitted by the secondary source
distributions, or even generate binaural room scanning (BRS) impulse response
sets in order to generate binaural simulations of the synthesized sound fields
with the `SoundScape Renderer`_.

.. _SoundScape Renderer: http://spatialaudio.net/ssr

Theory:
    http://sfstoolbox.org/

Documentation:
    http://matlab.sfstoolbox.org/

Source code and issue tracker:
    http://github.com/sfstoolbox/sfs-matlab/

SFS Toolbox for Python:
    http://python.sfstoolbox.org/

License:
    MIT -- see the file ``LICENSE`` for details.


Installation
------------

`Download the Toolbox`_, go to the main path of the Toolbox and start it with
``SFS_start`` which will add all needed paths to Matlab/Octave.  If
you want to remove them again, run ``SFS_stop``.

.. _Download the Toolbox: https://github.com/sfstoolbox/sfs-matlab/releases/latest


Requirements
------------

Matlab:
    You need Matlab version R2011b or newer to run the Toolbox.  On older
    versions the Toolbox should also work, but you need to add `narginchk.m`_ to
    the ``SFS_helper`` directory.

.. _narginchk.m: http://gist.github.com/hagenw/5642886

Octave:
    You need Octave version 3.6 or newer to run the Toolbox. In addition,
    you will need the ``audio`` and ``signal`` packages from
    `octave-forge`_.

.. _octave-forge: http://octave.sourceforge.net/

audioread:
    If ``audioread()`` is not available in your Matlab or Octave version,
    you can replace it by ``wavread()``. It is used in the two functions
    ``auralize_ir()`` and ``compensate_headphone()``.

Impulse responses:
    The Toolbox uses the `SOFA`_ file format for handling impulse response data
    sets like HRTFs. If you want to use this functionality you also have to
    install the `SOFA API for Matlab/Octave`_, which you can add to your paths
    by executing ``SOFAstart``.

Backward compatibility:
    Since version 2.0.0 the SFS Toolbox incorporates `SOFA`_ as file format for
    HRTFs which replaces the old `irs file format`_ formerly used by the
    Toolbox. If you still need this you should download `the latest version with
    irs file support`_.

.. _SOFA: http://sofaconventions.org/
.. _SOFA API for Matlab/Octave: https://github.com/sofacoustics/API_MO
.. _irs file format: https://dev.qu.tu-berlin.de/projects/measurements/wiki/IRs_file_format
.. _the latest version with irs file support: https://github.com/sfstoolbox/sfs-matlab/releases/tag/1.2.0


Getting started
---------------

In order to make a simulation of the sound field of a monochromatic point source
with a frequency of 800 Hz placed at (0,2.5,0) m synthesized by WFS run

.. sourcecode:: matlab

    conf = SFS_config;
    conf.plot.normalisation = 'center';
    sound_field_mono_wfs([-2 2],[-2 2],0,[0 2.5 0],'ps',800,conf)

To make a simulation of the same point source - now producing a broadband
impulse - in the time domain at a time of 5 ms after the first loudspeaker
activity run

.. sourcecode:: matlab

    conf = SFS_config;
    conf.plot.normalisation = 'max';
    sound_field_imp_wfs([-2 2],[-2 2],0,[0 2.5 0],'ps',0.005,conf)

After that have a look at ``SFS_config.m`` for the default settings of
the Toolbox.  Please don't change the settings directly in
``SFS_config.m``, but create an extra function or script for this, that
can look like this:

.. sourcecode:: matlab

    conf = SFS_config;
    conf.fs = 48000;

For a detailed description of all available features the SFS Toolbox, have a
look at the `online documentation`_.

.. _online documentation: http://matlab.sfstoolbox.org


Credits and feedback
--------------------

If you have questions, bug reports or feature requests, please use the `Issue
Section`_ to report them.

If you use the SFS Toolbox for your publications please cite our AES Convention
e-Brief and the DOI for the used Toolbox version, you will find at the `official
releases page`_:  

H. Wierstorf, S. Spors - Sound Field Synthesis Toolbox.
In the Proceedings of *132nd Convention of the
Audio Engineering Society*, 2012
[ `pdf`_ ]
[ `bibtex`_ ]

Copyright (c) 2010-2017 SFS Toolbox Developers

.. _Issue Section: https://github.com/sfstoolbox/sfs-matlab/issues
.. _official releases page: https://github.com/sfstoolbox/sfs-matlab/releases
.. _pdf: http://files.sfstoolbox.org/wierstorf_et_al_sfs-toolbox_aes132.pdf
.. _bibtex: http://files.sfstoolbox.org/wierstorf_et_al_sfs-toolbox_aes132.bib
