function status = test_binaural_synthesis(modus)
%TEST_BINAURAL_SYNTHESIS tests the correctness of the binaural synthesis
%functions
%
%   Usage: status = test_binaural_synthesis(modus)
%
%   Input parameters:
%       modus   - 0: numerical
%                 1: visual (not available)
%
%   Output parameters:
%       status  - true or false

%   TEST_BINAURAL_SYNTHESIS(modus) tests the ir_wfs function for different
%   loudspeaker arrays and source models.

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
%*****************************************************************************


status = false;


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
if modus
    warning('%s: visual modus not available.',upper(mfilename));
end


%% ===== Settings ========================================================
conf = SFS_config;
conf.c = 343;
conf.fs = 44100;
conf.secondary_sources.x0 = [];
conf.secondary_sources.center = [0 0 0];
conf.secondary_sources.size = 3;
conf.ir.usehcomp = false;
conf.ir.useinterpolation = true;
conf.N = 4096;
conf.dimension = '2.5D';
conf.driving_functions = 'default';
conf.usetapwin = true;
conf.tapwinlen = 0.3;
conf.wfs.usehpre = true;
conf.wfs.hpretype = 'FIR';
conf.wfs.hpreflow = 50;
conf.wfs.hprefhigh = 1200;
conf.delayline.resampling = 'none';
conf.delayline.filter = 'integer';
conf.debug = 0;
conf.showprogress = false;


%% ===== Main ============================================================
% Check if HRTF data set is available, download otherwise
basepath = get_sfs_path();
hrtf_file = [basepath '/data/HRTFs/QU_KEMAR_anechoic_3m.sofa'];
if ~exist(hrtf_file,'file')
    disp('Download');
    url = ['https://github.com/sfstoolbox/data/blob/master/', ...
           'HRTFs/QU_KEMAR_anechoic_3m.sofa?raw=true'];
    download_file(url,hrtf_file);
end
% Load HRTF data set
hrtf = SOFAload(hrtf_file);


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


status = true;
