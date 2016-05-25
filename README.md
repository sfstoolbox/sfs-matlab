Sound Field Synthesis Toolbox
=============================

Matlab/Octave implementation of the [Sound Field Synthesis
Toolbox](http://sfstoolbox.org).

The Sound Field Synthesis (SFS) Toolbox for Matlab/Octave gives you the
possibility to play around with sound field synthesis methods like Wave Field
Synthesis (WFS) or near-field compensated Higher Order Ambisonics (NFC-HOA).
There are functions to simulate monochromatic sound fields for different secondary
source (loudspeaker) setups, time snapshots of full band impulses emitted by the
secondary source distributions, or even generate Binaural Room Scanning (BRS)
stimuli sets in order to simulate WFS with the SoundScape Renderer (SSR).

**Documentation**
http://sfstoolbox.org/doc/latest/matlab/

**Mathematics**
http://sfstoolbox.org/doc/latest/math/

**License**
MIT -- see the file `LICENSE` for details.


Installation
------------

Download the Toolbox, go to the main path of the Toolbox and start it with
<code>SFS_start</code> which will add all needed paths to Matlab/Octave.  If
you want to remove them again, run <code>SFS_stop</code>.


Requirements
------------

**Matlab**
You need Matlab version 2011b or newer to run the Toolbox.
On older version the Toolbox should also work, but you need to add
[narginchk.m](http://gist.github.com/hagenw/5642886) to the
<code>SFS_helper</code>
directory.

**Octave**
You need Octave version 3.6 or newer to run the Toolbox. In addition,
you will need the following additional packages from
[octave-forge](http://octave.sourceforge.net/):
* audio
* signal (e.g. for firls)

**Impulse responses**
The Toolbox uses the [SOFA](http://sofaconventions.org/) file format for
handling impulse response data sets like HRTFs. If you want to use this
functionality you also have to install the [SOFA API for
Matlab/Octave](https://github.com/sofacoustics/API_MO), which you can add to
your paths by executing `SOFAstart`.

After setting up the Toolbox you can made one of the magic following things with it.


Getting started
---------------

In order to make a simulation of the sound field of a monochromatic point source
with a frequency of 800 Hz placed at (0,2.5,0) m synthesized by WFS run

```Matlab
conf = SFS_config;
conf.plot.normalisation = 'center';
sound_field_mono_wfs([-2 2],[-2 2],0,[0 2.5 0],'ps',800,conf)
```

To make a simulation of the same point source - now producing a broadband
impulse - in the time domain at a time of
200 samples after the first loudspeaker activity run

```Matlab
conf = SFS_config;
conf.plot.normalisation = 'max';
sound_field_imp_wfs([-2 2],[-2 2],0,[0 2.5 0],'ps',200,conf)
```

After that have a look at <code>SFS_config.m</code> for the default settings of
the Toolbox.  Please don't change the settings directly in
<code>SFS_config.m</code>, but create an extra function or script for this, that
can look like this:

```Matlab
conf = SFS_config;
conf.fs = 48000;
```

For a detailed description of all available features of the SFS Toolbox, have a
look at the [online documentation](http://sfstoolbox.org/doc/latest/matlab/) at
the [SFS Toolbox web page](http://sfstoolbox.org).


Credits and feedback
-------------------

If you have questions, bug reports or feature requests, please use the [Issue
Section](https://github.com/sfstoolbox/sfs/issues) to report them.

If you use the SFS Toolbox for your publications please cite our AES Convention
e-Brief and the DOI for the used Toolbox version, you will find at the [official
releases page](https://github.com/sfstoolbox/sfs/releases):  
H. Wierstorf, S. Spors - Sound Field Synthesis Toolbox.
In the Proceedings of *132nd Convention of the
Audio Engineering Society*, 2012
[ [pdf](http://audio.qu.tu-berlin.de/wp-content/uploads/publications/2012/wierstorf2012_SFS_toolbox_AES.pdf) ]
[ [bibtex](http://sfstoolbox.org/files/aes132_paper.bib) ]

Copyright (c) 2010-2016 SFS Toolbox Developers
