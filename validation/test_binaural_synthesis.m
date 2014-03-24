function test_binaural_synthesis()
%TEST_BINAURAL_SYNTHESIS tests the correctness of the binaural synthesis
%functions
%
%   Usage: test_binaural_synthesis()
%
%   TEST_BINAURAL_SYNTHESIS() tests the ir_wfs_25d function for different
%   loudspeaker arrays and source models.

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
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
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ===================================
nargmin = 0;
nargmax = 0;
narginchk(nargmin,nargmax);


%% ===== Settings ========================================================
conf.c = 343;
conf.fs = 44100;
conf.secondary_sources.x0 = [];
conf.secondary_sources.center = [0 0 0];
conf.secondary_sources.size = 3;
conf.ir.useoriglength = false;
conf.ir.usehcomp = false;
conf.ir.useinterpolation = true;
conf.N = 2048;
conf.dimension = '2.5D';
conf.driving_functions = 'default';
conf.usetapwin = true;
conf.tapwinlen = 0.3;
conf.wfs.usehpre = true;
conf.wfs.hpretype = 'FIR';
conf.wfs.hpreflow = 50;
conf.wfs.hprefhigh = 1200;
conf.usefracdelay = false;
conf.debug = 0;
conf.showprogress = false;


%% ===== Main ============================================================
% check if HRTF data set is available, download otherwise
basepath = get_sfs_path();
hrtf_file = [basepath '/data/HRTFs/QU_KEMAR_anechoic_3m.mat'];
if ~exist(hrtf_file,'file')
    url = ['https://dev.qu.tu-berlin.de/projects/measurements/repository/', ...
        'raw/2010-11-kemar-anechoic/mat/QU_KEMAR_anechoic_3m.mat'];
    download_file(url,hrtf_file);
end
% load HRTF data set
hrtf = read_irs(hrtf_file,conf);


%% ===== WFS 2.5D ========================================================
% === Linear array ===
conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.number = 20;
X = [0 -2 0];
phi = pi/2;
conf.xref = X;
% Plane wave
src = 'pw';
xs = [0.5 -1 0];
ir_wfs(X,phi,xs,src,hrtf,conf);
% Point source
src = 'ps';
xs = [0 1 0];
ir_wfs(X,phi,xs,src,hrtf,conf);
% Focused source
src = 'fs';
xs = [0 -1 0 0 -1 0];
ir_wfs(X,phi,xs,src,hrtf,conf);

% === Circular array ===
conf.secondary_sources.geometry = 'circle';
conf.secondary_sources.number = 64;
X = [0 0 0];
phi = pi/2;
conf.xref = X;
% Plane wave
src = 'pw';
xs = [0.5 -1 0];
ir_wfs(X,phi,xs,src,hrtf,conf);
% Point source
src = 'ps';
xs = [0.5 2 0];
ir_wfs(X,phi,xs,src,hrtf,conf);
% Focused source
src = 'fs';
xs = [0.5 0.5 0 -1 -1 0];
ir_wfs(X,phi,xs,src,hrtf,conf);

% === Box shaped array ===
conf.secondary_sources.geometry = 'box';
conf.secondary_sources.number = 80;
X = [0 0 0];
phi = pi/2;
conf.xref = X;
% Plane wave
src = 'pw';
xs = [0.5 -1 0];
ir_wfs(X,phi,xs,src,hrtf,conf);
% Point source
src = 'ps';
xs = [0.5 2 0];
ir_wfs(X,phi,xs,src,hrtf,conf);
% Focused source
src = 'fs';
xs = [0.5 0.5 0 -1 -1 0];
ir_wfs(X,phi,xs,src,hrtf,conf);
