function status = test_localwfs_vss(modus)
%TEST_LOCALWFS_VSS tests the LWFS-VSS driving functions
%
%   Usage: status = test_localwfs(modus)
%
%   Input parameters:
%       modus   - 0: numerical
%                 1: visual
%
%   Output parameters:
%       status  - true or false
%
%   TEST_LOCALWFS_VSS(modus) checks if the monochromatic and time-domain 
%   LWFS-VSS driving functions are working. Different sound fields are 
%   calculated and plotted for visual inspection.

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


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);


%% ===== Configuration ===================================================
% Parameters
conf = SFS_config;

conf.showprogress = true;

conf.plot.useplot = modus;
conf.plot.normalisation = 'center';
conf.plot.loudspeakers = true;
conf.plot.realloudspeakers = false;
conf.resolution = 500;
conf.t0 = 'source';

% default setting for WFS
conf.usetapwin = true;
conf.tapwinlen = 0.3;

% default settings for LWFS-VSS
conf.localwfs_vss.method = 'wfs';
conf.localwfs_vss.usetapwin = true;
conf.localsfstapwinlen = 0.3;
conf.localwfs_vss.consider_target_field = true;
conf.localwfs_vss.consider_secondary_sources = true;

% range for sound field computation
X = [-1 1];
Y = [-1 1];
Z = 0;

f = 7000;

conf.xref = [0,0,0];

% test scenarios
scenarios = { ...
    '2.5D', 'circular', 'circular', 'pw', [1 -1.0 0], 'mono'
    '2.5D', 'circular',   'linear', 'pw', [1 -1.0 0], 'mono'
    '2.5D', 'circular',      'box', 'pw', [1 -1.0 0], 'mono'
    '2.5D',   'linear', 'circular', 'pw', [1 -1.0 0], 'mono'
    '2.5D',   'linear',   'linear', 'pw', [1 -1.0 0], 'mono'
    '2.5D',   'linear',      'box', 'pw', [1 -1.0 0], 'mono'
    '2.5D',      'box', 'circular', 'pw', [1 -1.0 0], 'mono'
    '2.5D',      'box',   'linear', 'pw', [1 -1.0 0], 'mono'
    '2.5D',      'box',      'box', 'pw', [1 -1.0 0], 'mono'
    '2.5D', 'circular', 'circular', 'ps',  [1 2.5 0], 'mono'
    '2.5D', 'circular',   'linear', 'ps',  [1 2.5 0], 'mono'
    '2.5D', 'circular',      'box', 'ps',  [1 2.5 0], 'mono'
    '2.5D',   'linear', 'circular', 'ps',  [1 2.5 0], 'mono'
    '2.5D',   'linear',   'linear', 'ps',  [1 2.5 0], 'mono'
    '2.5D',   'linear',      'box', 'ps',  [1 2.5 0], 'mono'
    '2.5D',      'box', 'circular', 'ps',  [1 2.5 0], 'mono'
    '2.5D',      'box',   'linear', 'ps',  [1 2.5 0], 'mono'
    '2.5D',      'box',      'box', 'ps',  [1 2.5 0], 'mono'
    '2.5D', 'circular', 'circular', 'pw', [1 -1.0 0], 'imp'
    '2.5D', 'circular',   'linear', 'pw', [1 -1.0 0], 'imp'
    '2.5D', 'circular',      'box', 'pw', [1 -1.0 0], 'imp'
    '2.5D',   'linear', 'circular', 'pw', [1 -1.0 0], 'imp'
    '2.5D',   'linear',   'linear', 'pw', [1 -1.0 0], 'imp'
    '2.5D',   'linear',      'box', 'pw', [1 -1.0 0], 'imp'
    '2.5D',      'box', 'circular', 'pw', [1 -1.0 0], 'imp'
    '2.5D',      'box',   'linear', 'pw', [1 -1.0 0], 'imp'
    '2.5D',      'box',      'box', 'pw', [1 -1.0 0], 'imp'
    '2.5D', 'circular', 'circular', 'ps',  [1 2.5 0], 'imp'
    '2.5D', 'circular',   'linear', 'ps',  [1 2.5 0], 'imp'
    '2.5D', 'circular',      'box', 'ps',  [1 2.5 0], 'imp'
    '2.5D',   'linear', 'circular', 'ps',  [1 2.5 0], 'imp'
    '2.5D',   'linear',   'linear', 'ps',  [1 2.5 0], 'imp'
    '2.5D',   'linear',      'box', 'ps',  [1 2.5 0], 'imp'
    '2.5D',      'box', 'circular', 'ps',  [1 2.5 0], 'imp'
    '2.5D',      'box',   'linear', 'ps',  [1 2.5 0], 'imp'
    '2.5D',      'box',      'box', 'ps',  [1 2.5 0], 'imp'
    };

% Start testing
for ii=1:size(scenarios)
    
    conf.dimension = scenarios{ii,1};
    
    % set loudspeaker array
    switch scenarios{ii,2}
    case 'linear'
        conf.secondary_sources.geometry = 'linear';
        conf.secondary_sources.number = 56;
        conf.secondary_sources.size = 2;
        conf.secondary_sources.center = [0 1 0];
    case 'circular'
        conf.secondary_sources.geometry = 'circular';
        conf.secondary_sources.number = 56;
        conf.secondary_sources.size = 2;
        conf.secondary_sources.center = [0 0 0];
    case 'box'
        conf.secondary_sources.geometry = 'box';
        conf.secondary_sources.number = 56*4;
        conf.secondary_sources.size = 2;
        conf.secondary_sources.center = [0 0 0];
    end
    
    % set virtual source distribution
    switch scenarios{ii,3}
    case 'linear'
        conf.localwfs_vss.size = 0.4;
        conf.localwfs_vss.center = [0 0.2 0];
        conf.localwfs_vss.geometry = 'linear';
        conf.localwfs_vss.number = 56;
    case 'circular'
        conf.localwfs_vss.size = 0.4;
        conf.localwfs_vss.center = [0 0 0];
        conf.localwfs_vss.geometry = 'circular';
        conf.localwfs_vss.number = 56;
    case 'box'
        conf.localwfs_vss.size = 0.4;
        conf.localwfs_vss.center = [0 0 0];
        conf.localwfs_vss.geometry = 'box';
        conf.localwfs_vss.number = 4*56;
    end
    
    % set desired sound field
    src = scenarios{ii,4};
    xs = scenarios{ii,5};
    
    % set domain (monochromatic or time)
    switch scenarios{ii,6}
    case 'mono'
        conf.plot.usedb = false;
        [~] = sound_field_mono_localwfs_vss(X,Y,Z,xs,src,f,conf);
    case 'imp'
        % set t for time-shapshot
        switch src
        case {'ps', 'ls', 'fs'}
            t = norm(xs - conf.xref)/conf.c;
        case {'pw'}
            t = 0;
        otherwise
            error('unknown source type');
        end
        conf.plot.usedb = true;
        [~] = sound_field_imp_localwfs_vss(X,Y,Z,xs,src,t,conf);
    end
    
    % title of plot
    if modus
        title(sprintf(['LWFS-VSS (%s,%s): %s SSD, %s virtual SSD, %s' ...
            ' (%02.02f, %02.02f, %02.02f)'], scenarios{ii,1}, scenarios{ii,6}...
            , scenarios{ii,2}, scenarios{ii,3}, src, xs));
    end
end

status = true;
