function fal = aliasing_frequency(dx0,conf)
%ALIASING_FREQUENCY returns the aliasing frequency
%
%   Usage: fal = aliasing_frequency(dx0,[conf])
%
%   Input options:
%       dx0     - distance between adjacent secondary sources (m)
%       conf    - optional configuration struct (see SFS_config)
%
%   Output options:
%       fal     - aliasing frequency (Hz)
%
%   ALIASING_FREQUENCY(dx0,conf) returns the aliasing frequency for the given
%   interspacing of secondary sources. The value is calculated after
%   spors2009.
%   
%   References:
%       S. Spors and J. Ahrens - Spatial sampling artifacts of wave field
%       synthesis for the reproduction of virtual point sources. 126th AES,
%       May 2009.
%       E. Start - Direct Sound Enhancement by Wave Field Synthesis. TU Delft,
%       1997.
%
%   see also: wave_field_mono_wfs_25d

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
if nargin==nargmax-1
    conf = SFS_config;
end
%[xs,X] = position_vector(xs,X);


%% ===== Configuration ==================================================
c = conf.c;


%% ===== Computation =====================================================
fal = c/(2*dx0);
return

%FIXME: the following code works not probably at the moment.
% The alias frequency depends on the plane waves that are present in the sound
% field. And which plane waves can be present depends on the concrete angle
% between the source/listener and the array of secondary sources. The following
% code uses the Fraunhofer approximation in order to calculate the alias
% frequency, see Start1997 p73, eq 5.17
%
% get secondary source positions and activity
x0 = secondary_source_positions(L,conf);
x0 = secondary_source_selection(x0,xs,src);
win = tapering_window(x0,conf);
% === get angles between array edges and virtual source position ===
% orientation of loudspeaker array
v_array = x0(1,1:3)-x0(end,1:3);
[phi_array,theta_tmp,r_tmp] = cart2sph(v_array(1),v_array(2),v_array(3))-pi/2;
% angle between array edges and virtual source position
v1 = x0(1,1:3)-xs;
v2 = x0(end,1:3)-xs;
[phix01,theta_tmp,r_tmp] = cart2sph(v1(1),v1(2),v1(3))-phi_array;
[phix02,theta_tmp,r_tmp] = cart2sph(v2(1),v2(2),v2(3))-phi_array;
% angle between array edges and listener position
v1 = x0(1,1:3)-X;
v2 = x0(end,1:3)-X;
[phiX1,theta_tmp,r_tmp] = cart2sph(v1(1),v1(2),v1(3))+phi_array;
[phiX2,theta_tmp,r_tmp] = cart2sph(v2(1),v2(2),v2(3))+phi_array;
% use only the smallest of these angles
phi1 = min([abs(phix01) abs(phiX1)]);
phi2 = min([abs(phix02) abs(phiX2)]);
%                    c
% fal = -----------------------------
%       dx0 (sin(|phi1|)+sin(|phi2|))
fal = c / (dx0*(sin(phi1)+sin(phi2)));
