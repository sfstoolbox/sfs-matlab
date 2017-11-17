.. _sec-usage:

Usage
=====


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


Installation
------------

`Download the Toolbox`_, go to the main path of the Toolbox and start it with
``SFS_start`` which will add all needed paths to Matlab/Octave.  If
you want to remove them again, run ``SFS_stop``.

.. _Download the Toolbox: https://github.com/sfstoolbox/sfs-matlab/releases/latest


How to Get Started
------------------

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
