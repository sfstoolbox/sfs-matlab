Sound Field Synthesis Toolbox
=============================

The Sound Field Synthesis Toolbox (SFS) for Matlab/Octave gives you the
possibility to play around with sound field synthesis methods like Wave Field
Synthesis (WFS) or near-field compensated Higher Order Ambisonics (NFCHOA).
There are functions to simulate monochromatic wave fields for different secondary
source (loudspeaker) setups, time snapshots of full band impulses emitted by the
secondary source distributions, or even generate Binaural Room Scanning (BRS)
stimuli sets in order to simulate WFS with the SoundScape Renderer (SSR).


Installation
------------

Download the Toolbox and add the main path of the Toolbox to your Matlab/Octave
path. After that copy <code>SFS_config_example.m</code> to
<code>SFS_config.m</code> and change it to your needs. For an easy beginning you
can just use the default settings by leaving everything as it is.
Then start the Toolbox with <code>SFS_start</code> which will add all needed
subpathes.


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
* audio (e.g. for wavwrite)
* signal (e.g. for firls)

Now you set up the Toolbox and can made on of the following things with it:

Secondary Sources
-----------------

The Toolbox comes with a function which can generate different common shapes of loudspeaker arrays for you.
At the moment these include linear, circular and box shaped arrays.

Before showing the different geometries, we start with some common settings. First we get a configuration struct
and set the array size/diameter to 3m.

```Matlab
conf = SFS_config;
L = 3;
```

### linear array

```Matlab
```

![Image](doc/img/secondary_sources_linear.png)

Simulate monochromatic sound fields
-----------------------------------

With the files in <code>SFS_monochromatic</code> you can simulate a
monochromatic sound field in a specified area for different techniques like WFS
and NFCHOA.


Simulate time snapshots of sound fields
---------------------------------------

With the files in <code>SFS_time_domain</code> you can simulate snapshots in
time of an impulse sendiong out from your WFS or NFCHOA system


Make binaural simulations of your systems
-----------------------------------------

If you have a set of head-related transfer functions (HRTFs) you can simulate
the ear signals reaching a listener sitting at a given point in the listening
area for a specified WFS or NFCHOA system.
You can even download a set of HRTFs, which will just work with the Toolbox at 
http://dev.qu.tu-berlin.de/projects/measurements/wiki/2010-11-kemar-anechoic

The files dealing with the binaural simulations are in the folder
<code>SFS_binaural_synthesis</code>. Files dealing with HRTFs are in the folder
<code>SFS_ir</code>. If you want to extrapolate your HRTFs to plane waves you
may also want to have a look in <code>SFS_HRTF_extrapolation</code>.


Small helper functions
----------------------

The Toolbox provides you also with a set of useful small functions that may want
to use. Here the highlights are angle conversion with <code>rad()</code> and
<code>degree()</code>, FFT calculation and plotting <code>easyfft()</code>,
create noise signal <code>noise()</code>, rotation matrix
<code>rotation_matrix()</code>, even or odd checking <code>iseven()</code>
<code>isodd()</code>, spherical bessel functions <code>sphbesselh()</code>
<code>sphbesselj</code> <code>sphbessely</code>.


Plotting with Matlab or Gnuplot
-------------------------------

The Toolbox provides you with a variety of functions for plotting your simulated
sound fields <code>plot_wavefield()</code> and adding loudspeaker symbols to the
figure <code>draw_loudspeakers</code>. If you have gnuplot installed, you can
even use it with the Toolbox by setting <code>conf.plot.usegnuplot =
true;</code>.


Credits and License
-------------------

This is the source distribution of Sound Field Synthesis Toolbox (SFS) licensed
under the GPLv3+. Please consult the file COPYING for more information about
this license.
 
For questions, bug reports and feature requests:  
Contact: sfstoolbox@googlemail.com  
Website: http://github.com/sfstoolbox/sfs


Copyright (c) 2010-2013  
Quality & Usability Lab, together with  
Assessment of IP-based Applications  
Telekom Innovation Laboratories, TU Berlin  
Ernst-Reuter-Platz 7, 10587 Berlin, Germany 


Copyright (c) 2013  
Institut fuer Nachrichtentechnik  
Universitaet Rostock  
Richard-Wagner-Strasse 31, 18119 Rostock
