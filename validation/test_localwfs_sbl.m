function status = test_localwfs_sbl(modus)
%TEST_LOCALWFS_SBL tests the LWFS-SBL driving functions
%
%   Usage: status = test_localwfs_sbl(modus)
%
%   Input parameters:
%       modus   - 0: numerical
%                 1: visual
%
%   Output parameters:
%       status  - true or false

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2019 SFS Toolbox Developers                             *
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
% https://sfs.readthedocs.io                            sfstoolbox@gmail.com *
%*****************************************************************************


status = false;


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);


%% ===== Configuration ===================================================
% Parameters
conf = SFS_config;

conf.plot.useplot = modus;
conf.plot.normalisation = 'center';
conf.plot.loudspeakers = true;
conf.plot.realloudspeakers = false;
conf.resolution = 500;
conf.t0 = 'source';

% Default setting for WFS
conf.usetapwin = true;
conf.tapwinlen = 0.3;
% Default settings for LWFS-SBL
conf.localwfs_sbl.fc = 1200;
conf.localwfs_sbl.Npw = 1024;

% Range for sound field computation
X = [-1 1];
Y = [-1 1];
Z = 0;

f = 4000;

conf.xref = [0,0,0];

% Test scenarios
scenarios = { ...
    '2.5D', 'circular', 10, 'pw', [1 -1.0 0], 'mono'
    '2.5D', 'circular', 20, 'pw', [1 -1.0 0], 'mono'
    '2.5D', 'circular', 30, 'pw', [1 -1.0 0], 'mono'
    '2.5D',   'linear', 10, 'pw', [1 -1.0 0], 'mono'
    '2.5D',   'linear', 20, 'pw', [1 -1.0 0], 'mono'
    '2.5D',   'linear', 30, 'pw', [1 -1.0 0], 'mono'
    '2.5D',      'box', 10, 'pw', [1 -1.0 0], 'mono'
    '2.5D',      'box', 20, 'pw', [1 -1.0 0], 'mono'
    '2.5D',      'box', 30, 'pw', [1 -1.0 0], 'mono'
    '2.5D', 'circular', 10, 'ps',  [1 2.5 0], 'mono'
    '2.5D', 'circular', 20, 'ps',  [1 2.5 0], 'mono'
    '2.5D', 'circular', 30, 'ps',  [1 2.5 0], 'mono'
    '2.5D',   'linear', 10, 'ps',  [1 2.5 0], 'mono'
    '2.5D',   'linear', 20, 'ps',  [1 2.5 0], 'mono'
    '2.5D',   'linear', 30, 'ps',  [1 2.5 0], 'mono'
    '2.5D',      'box', 10, 'ps',  [1 2.5 0], 'mono'
    '2.5D',      'box', 20, 'ps',  [1 2.5 0], 'mono'
    '2.5D',      'box', 30, 'ps',  [1 2.5 0], 'mono'
    '2.5D', 'circular', 10, 'pw', [1 -1.0 0], 'imp'
    '2.5D', 'circular', 20, 'pw', [1 -1.0 0], 'imp'
    '2.5D', 'circular', 30, 'pw', [1 -1.0 0], 'imp'
    '2.5D',   'linear', 10, 'pw', [1 -1.0 0], 'imp'
    '2.5D',   'linear', 20, 'pw', [1 -1.0 0], 'imp'
    '2.5D',   'linear', 30, 'pw', [1 -1.0 0], 'imp'
    '2.5D',      'box', 10, 'pw', [1 -1.0 0], 'imp'
    '2.5D',      'box', 20, 'pw', [1 -1.0 0], 'imp'
    '2.5D',      'box', 30, 'pw', [1 -1.0 0], 'imp'
    '2.5D', 'circular', 10, 'ps',  [1 2.5 0], 'imp'
    '2.5D', 'circular', 20, 'ps',  [1 2.5 0], 'imp'
    '2.5D', 'circular', 30, 'ps',  [1 2.5 0], 'imp'
    '2.5D',   'linear', 10, 'ps',  [1 2.5 0], 'imp'
    '2.5D',   'linear', 20, 'ps',  [1 2.5 0], 'imp'
    '2.5D',   'linear', 30, 'ps',  [1 2.5 0], 'imp'
    '2.5D',      'box', 10, 'ps',  [1 2.5 0], 'imp'
    '2.5D',      'box', 20, 'ps',  [1 2.5 0], 'imp'
    '2.5D',      'box', 30, 'ps',  [1 2.5 0], 'imp'
    };

% Start testing
for ii=1:size(scenarios,1)
    
    conf.dimension = scenarios{ii,1};
    
    % Set loudspeaker array
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
    
    conf.localwfs_sbl.order = scenarios{ii,3};
    
    % Set desired sound field
    src = scenarios{ii,4};
    xs = scenarios{ii,5};
    
    % Set domain (monochromatic or time)
    switch scenarios{ii,6}
    case 'mono'
        conf.plot.usedb = false;
        [~] = sound_field_mono_localwfs_sbl(X,Y,Z,xs,src,f,conf);
    case 'imp'
        % Set t for time-shapshot
        switch src
        case {'ps', 'ls', 'fs'}
            t = norm(xs - conf.xref)/conf.c;
        case {'pw'}
            t = 0;
        otherwise
            error('unknown source type');
        end
        conf.plot.usedb = true;
        [~] = sound_field_imp_localwfs_sbl(X,Y,Z,xs,src,t,conf);
    end
    
    progress_bar(ii,size(scenarios,1));

    % Title of plot
    if modus
        title(sprintf(['LWFS-SBL (%s,%s): %s SSD, order: %d, %s' ...
            ' (%02.02f, %02.02f, %02.02f)'], scenarios{ii,1}, scenarios{ii,6}...
            , scenarios{ii,2}, scenarios{ii,3}, src, xs));
    end
end


status = true;
