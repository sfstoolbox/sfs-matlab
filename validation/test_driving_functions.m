function boolean = test_driving_functions(modus)
%TEST_DRIVING_FUNCTIONS tests the correctness of the driving functions
%
%   Usage: boolean = test_driving_functions(modus)
%
%   Input parameters:
%       modus   - 0: numerical
%                 1: visual
%
%   Output parameters:
%       booelan - true or false
%
%   TEST_DRIVING_FUNCTIONS(MODUS) checks if the functions, that calculates
%   the driving functions working correctly. Therefore different wave
%   fields are simulated.

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

% AUTHOR: Hagen Wierstorf
% $LastChangedDate: $
% $LastChangedRevision: $
% $LastChangedBy: $

% TODO: add mode to save data as reference data


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));


%% ===== Main ============================================================
conf = SFS_config;
L = 3;
f = 1000;
conf.useplot = modus;


%% ===== WFS 2D ==========================================================
% === Linear array ===
conf.array = 'linear';
conf.dx0 = 0.15;
X = [-2,2];
Y = [-0.15,3];
conf.xref = [0 2];
% Plane wave
src = 'pw';
xs = [0.5,1];
% mono-frequent
[x,y,P_wfs2d_linear_pw,ls_activity] = wave_field_mono_wfs_2d(X,Y,xs,src,f,L,conf);
title('WFS 2D linear array, plane wave, mono-frequent');
% spatio-temporal impulse response (not implemented yet)
%[x,y,p_wfs2d_linear_pw,ls_activity] = wave_field_imp_wfs_2d(X,Y,xs,src,L,conf);
%title('WFS 2D linear array, plane wave, impulse response');
% Line source
src = 'ls';
xs = [0,-1];
% mono-frequent
[x,y,P_wfs2d_linear_ps,ls_activity] = wave_field_mono_wfs_2d(X,Y,xs,src,f,L,conf);
title('WFS 2D linear array, line source, mono-frequent');
% spatio-temporal impulse response (not implemented yet)
%[x,y,p_wfs2d_linear_ps,ls_activity] = wave_field_imp_wfs_2d(X,Y,xs,src,L,conf);
%title('WFS 2D linear array, line source, impulse response');
% Focused source
src = 'fs';
xs = [0,1];
% mono-frequent
[x,y,P_wfs2d_linear_fs,ls_activity] = wave_field_mono_wfs_2d(X,Y,xs,src,f,L,conf);
title('WFS 2D linear array, focused source, mono-frequent');
% spatio-temporal impulse response (not implemented yet)
%[x,y,p_wfs2d_linear_fs,ls_activity] = wave_field_imp_wfs_2d(X,Y,xs,src,L,conf);
%title('WFS 2D linear array, focused source, impulse response');

% === Circular array ===
conf.array = 'circle';
conf.dx0 = 0.15;
X = [-2,2];
Y = [-2,2];
conf.xref = [0,0];
% Plane wave
src = 'pw';
xs = [0.5,1];
% mono-frequent
[x,y,P_wfs2d_circular_pw,ls_activity] = wave_field_mono_wfs_2d(X,Y,xs,src,f,L,conf);
title('WFS 2D circular array, plane wave, mono-frequent');
% spatio-temporal impulse response (not implemented yet)
%[x,y,p_wfs2d_circular_pw,ls_activity] = wave_field_imp_wfs_2d(X,Y,xs,src,L,conf);
%title('WFS 2D circular array, plane wave, impulse response');
% Line source
src = 'ls';
xs = [0.5,2];
% mono-frequent
[x,y,P_wfs2d_circular_ps,ls_activity] = wave_field_mono_wfs_2d(X,Y,xs,src,f,L,conf);
title('WFS 2D circular array, line source, mono-frequent');
% spatio-temporal impulse response (not implemented yet)
%[x,y,p_wfs2d_circular_ps,ls_activity] = wave_field_imp_wfs_2d(X,Y,xs,src,L,conf);
%title('WFS 2D circular array, line source, impulse response');
% Focused source
src = 'fs';
xs = [0.5,0.5];
% mono-frequent
[x,y,P_wfs2d_circular_fs,ls_activity] = wave_field_mono_wfs_2d(X,Y,xs,src,f,L,conf);
title('WFS 2D circular array, focused source, mono-frequent');
% spatio-temporal impulse response (not implemented yet)
%[x,y,p_wfs2d_circular_fs,ls_activity] = wave_field_imp_wfs_2d(X,Y,xs,src,L,conf);
%title('WFS 2D circular array, focused source, impulse response');

% === Box shaped array ===
conf.array = 'box';
conf.dx0 = 0.15;
X = [-2,2];
Y = [-2,2];
conf.xref = [0,0];
% Plane wave
src = 'pw';
xs = [0.5,1];
% mono-frequent
[x,y,P_wfs2d_box_pw,ls_activity] = wave_field_mono_wfs_2d(X,Y,xs,src,f,L,conf);
title('WFS 2D box shaped array, plane wave, mono-frequent');
% spatio-temporal impulse response (not implemented yet)
%[x,y,p_wfs2d_box_pw,ls_activity] = wave_field_imp_wfs_2d(X,Y,xs,src,L,conf);
%title('WFS 2D box shaped array, plane wave, impulse response');
% Line source
src = 'ls';
xs = [0.5,2];
% mono-frequent
[x,y,P_wfs2d_box_ps,ls_activity] = wave_field_mono_wfs_2d(X,Y,xs,src,f,L,conf);
title('WFS 2D box shaped array, line source, mono-frequent');
% spatio-temporal impulse response (not implemented yet)
%[x,y,p_wfs2d_box_ps,ls_activity] = wave_field_imp_wfs_2d(X,Y,xs,src,L,conf);
%title('WFS 2D box shaped array, line source, impulse response');
% Focused source
src = 'fs';
xs = [0.5,0.5];
% mono-frequent
[x,y,P_wfs2d_box_fs,ls_activity] = wave_field_mono_wfs_2d(X,Y,xs,src,f,L,conf);
title('WFS 2D box shaped array, focused source, mono-frequent');
% spatio-temporal impulse response (not implemented yet)
%[x,y,p_wfs2d_box_fs,ls_activity] = wave_field_imp_wfs_2d(X,Y,xs,src,L,conf);
%title('WFS 2D box shaped array, focused source, impulse response');


%% ===== WFS 2.5D ========================================================
% === Linear array ===
conf.array = 'linear';
conf.dx0 = 0.15;
X = [-2,2];
Y = [-0.15,3];
conf.xref = [0 2];
conf.frame = 200;
% Plane wave
src = 'pw';
xs = [0.5,1];
% mono-frequent
[x,y,P_wfs25d_linear_pw,ls_activity] = wave_field_mono_wfs_25d(X,Y,xs,src,f,L,conf);
title('WFS 2.5D linear array, plane wave [0.5,1], mono-frequent');
% spatio-temporal impulse response
[x,y,p_wfs25d_linear_pw,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,src,L,conf);
title('WFS 2.5D linear array, plane wave [0.5,1], impulse response');
% Point source
src = 'ps';
xs = [0,-1];
% mono-frequent
[x,y,P_wfs25d_linear_ps,ls_activity] = wave_field_mono_wfs_25d(X,Y,xs,src,f,L,conf);
title('WFS 2.5D linear array, point source [0,-1], mono-frequent');
% spatio-temporal impulse response
[x,y,p_wfs25d_linear_ps,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,src,L,conf);
title('WFS 2.5D linear array, point source [0,-1], impulse response');
% Focused source
src = 'fs';
xs = [0,1];
% mono-frequent
[x,y,P_wfs25d_linear_fs,ls_activity] = wave_field_mono_wfs_25d(X,Y,xs,src,f,L,conf);
title('WFS 2.5D linear array, focused source [0,1], mono-frequent');
% spatio-temporal impulse response
[x,y,p_wfs25d_linear_fs,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,src,L,conf);
title('WFS 2.5D linear array, focused source [0,1], impulse response');

% === Circular array ===
conf.array = 'circle';
conf.dx0 = 0.15;
X = [-2,2];
Y = [-2,2];
conf.xref = [0,0];
% Plane wave
src = 'pw';
xs = [0.5,1];
% mono-frequent
[x,y,P_wfs25d_circular_pw,ls_activity] = wave_field_mono_wfs_25d(X,Y,xs,src,f,L,conf);
title('WFS 2.5D circular array, plane wave [0.5,1], mono-frequent');
% spatio-temporal impulse response
[x,y,p_wfs25d_circular_pw,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,src,L,conf);
title('WFS 2.5D circular array, plane wave [0.5,1], impulse response');
% Line source
src = 'ps';
xs = [0.5,2];
% mono-frequent
[x,y,P_wfs25d_circular_ps,ls_activity] = wave_field_mono_wfs_25d(X,Y,xs,src,f,L,conf);
title('WFS 2.5D circular array, point source [0.5,2], mono-frequent');
% spatio-temporal impulse response
[x,y,p_wfs25d_circular_ps,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,src,L,conf);
title('WFS 2.5D circular array, point source [0.5,2], impulse response');
% Focused source
src = 'fs';
xs = [0.5,0.5];
% mono-frequent
[x,y,P_wfs25d_circular_fs,ls_activity] = wave_field_mono_wfs_25d(X,Y,xs,src,f,L,conf);
title('WFS 2.5D circular array, focused source [0.5,0.5], mono-frequent');
% spatio-temporal impulse response
[x,y,p_wfs25d_circular_fs,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,src,L,conf);
title('WFS 2.5D circular array, focused source [0.5,0.5], impulse response');

% === Box shaped array ===
conf.array = 'box';
conf.dx0 = 0.15;
X = [-2,2];
Y = [-2,2];
conf.xref = [0,0];
% Plane wave
src = 'pw';
xs = [0.5,1];
% mono-frequent
[x,y,P_wfs25d_box_pw,ls_activity] = wave_field_mono_wfs_25d(X,Y,xs,src,f,L,conf);
title('WFS 2.5D box shaped array, plane wave [0.5,1], mono-frequent');
% spatio-temporal impulse response
[x,y,p_wfs25d_box_pw,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,src,L,conf);
title('WFS 2.5D box shaped array, plane wave [0.5,1], impulse response');
% point source
src = 'ps';
xs = [0.5,2];
% mono-frequent
[x,y,P_wfs25d_box_ps,ls_activity] = wave_field_mono_wfs_25d(X,Y,xs,src,f,L,conf);
title('WFS 2.5D box shaped array, point source [0.5,2], mono-frequent');
% spatio-temporal impulse response
[x,y,p_wfs25d_box_ps,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,src,L,conf);
title('WFS 2.5D box shaped array, point source [0.5,2], impulse response');
% Focused source
src = 'fs';
xs = [0.5,0.5];
% mono-frequent
[x,y,P_wfs25d_box_fs,ls_activity] = wave_field_mono_wfs_25d(X,Y,xs,src,f,L,conf);
title('WFS 2.5D box shaped array, focused source [0.5,0.5], mono-frequent');
% spatio-temporal impulse response
[x,y,p_wfs25d_box_fs,ls_activity] = wave_field_imp_wfs_25d(X,Y,xs,src,L,conf);
title('WFS 2.5D box shaped array, focused source [0.5,0.5], impulse response');


%% ===== NFC-HOA 2.5D ========================================================
% === Circular array ===
conf.array = 'circle';
conf.dx0 = 0.15;
X = [-2,2];
Y = [-2,2];
conf.xref = [0,0];
conf.frame = 200;
% Plane wave
src = 'pw';
xs = [0.5,1];
% mono-frequent
[x,y,P_nfchoa25d_circular_pw,ls_activity] = ...
    wave_field_mono_nfchoa_25d(X,Y,xs,src,f,L,conf);
title('NFC-HOA 2.5D circular array, plane wave [0.5,1], mono-frequent');
% spatio-temporal impulse response
[x,y,p_nfchoa25d_circular_pw] = wave_field_imp_nfchoa_25d(X,Y,xs,src,L,conf);
title('NFC-HOA 2.5D circular array, plane wave [0.5,1], impulse response');

% point source
src = 'ps';
xs = [0.5,2];
% mono-frequent
%[x,y,P_nfchoa25d_circular_ps,ls_activity] = ...
%    wave_field_mono_nfchoa_25d(X,Y,xs,src,f,L,conf);
%title('NFC-HOA 2.5D circular array, point source [0.5,2], mono-frequent');
% spatio-temporal impulse response
[x,y,p_nfchoa25d_circular_ps] = wave_field_imp_nfchoa_25d(X,Y,xs,src,L,conf);
title('NFC-HOA 2.5D circular array, point source [0.5,2], impulse response');


