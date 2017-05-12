function status = test_driving_functions_imp_with_delay(modus)
%TEST_DRIVING_FUNCTIONS_IMP_WITH_DELAY tests the correctness of the time-domain
%driving functions, in the case of conf.t0 = 'source'
%
%   Usage: status = test_driving_functions_imp_with_delay(modus)
%
%   Input parameters:
%       modus   - 0: numerical
%                 1: visual
%
%   Output parameters:
%       status - true or false
%
%   TEST_DRIVING_FUNCTIONS_IMP_WITH_DELAY(MODUS) checks if the functions,
%   that calculates the WFS driving functions in time-domain work correctly
%   for the setting conf.t0 = 'source'.
%   Therefore different sound fields are simulated.

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


% TODO: add mode to save data as reference data
status = false;


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);


%% ===== Configuration ===================================================
conf = SFS_config;
conf.t0 = 'source';
conf.secondary_sources.size = 3;
conf.plot.useplot = false;
conf.plot.usenormalisation = true;
conf.plot.usedb = true;
conf.driving_functions = 'default';
% test scenarios
scenarios = { ...
    'WFS', '2D',   'linear',   'pw', [ 0.5 -1.0  0.0]; ...
    'WFS', '2D',   'linear',   'ls', [ 0.0  1.0  0.0]; ...
    'WFS', '2D',   'linear',   'fs', [ 0.0 -1.0  0.0  0.0 -1.0  0.0]; ...
    'WFS', '2D',   'circular', 'pw', [ 0.5  0.5  0.0]; ...
    'WFS', '2D',   'circular', 'ls', [ 0.0  2.5  0.0]; ...
    'WFS', '2D',   'circular', 'ps', [ 0.0  2.5  0.0]; ...
    'WFS', '2D',   'circular', 'fs', [ 0.0  0.5  0.0  0.0 -1.0  0.0]; ...
    'WFS', '2D',   'box',      'pw', [ 0.5  1.0  0.0]; ...
    'WFS', '2D',   'box',      'ls', [ 2.0  2.0  0.0]; ...
    'WFS', '2D',   'box',      'ps', [ 2.0  2.0  0.0]; ...
    'WFS', '2D',   'box',      'fs', [ 0.5  0.5  0.0  -1.0  -1.0  0.0]; ...
    'WFS', '2.5D', 'linear',   'pw', [ 0.5 -1.0  0.0]; ...
    'WFS', '2.5D', 'linear',   'ps', [ 0.0  1.0  0.0]; ...
    'WFS', '2.5D', 'linear',   'fs', [ 0.0 -1.0  0.0  0.0 -1.0  0.0]; ...
    'WFS', '2.5D', 'circular', 'pw', [ 0.5  0.5  0.0]; ...
    'WFS', '2.5D', 'circular', 'ps', [ 0.0  2.5  0.0]; ...
    'WFS', '2.5D', 'circular', 'fs', [ 0.0  0.5  0.0  0.0 -1.0  0.0]; ...
    'WFS', '2.5D', 'box',      'pw', [ 0.5  1.0  0.0]; ...
    'WFS', '2.5D', 'box',      'ps', [ 2.0  2.0  0.0]; ...
    'WFS', '2.5D', 'box',      'fs', [ 0.5  0.5  0.0 -1.0  -1.0  0.0]; ...
    'WFS', '3D',   'sphere',   'ps', [ 0.0  0.0  2.5]; ...
    'WFS', '3D',   'sphere',   'ls', [ 0.0  0.0  2.5  0.0  1.0  1.0]; ...
    'WFS', '3D',   'sphere',   'pw', [ 0.0  0.0  -1.0]; ...
    'WFS', '3D',   'sphere',   'fs', [ 0.0  0.0  0.5  0.0  0.0  -1.0]; ...
    'HOA', '2.5D', 'circular', 'pw', [ 0.5  0.5  0.0]; ...
    'HOA',   '3D',   'sphere', 'pw', [ 0.0  0.0 -1.0]; ...
    'HOA',   '3D',   'sphere', 'pw', [-1.0  0.0  0.0]; ...
    'HOA', '2.5D', 'circular', 'ps', [ 0.0  2.5  0.0]; ...
    'HOA',   '3D',   'sphere', 'ps', [ 0.0  0.0  2.5]; ...
    'HOA',   '3D',   'sphere', 'ps', [ 2.5  0.0  0.0]; ...
};

% Start testing
for ii=1:size(scenarios)

    % get current dimension
    conf.dimension = scenarios{ii,2};

    % get listening area
    switch scenarios{ii,3}
        case 'linear'
            X = [-2 2];
            Y = [-3 0.15];
            Z = 0;
            conf.xref = [0 -1.5 0];
            conf.secondary_sources.number = 20;
        case 'circular'
            X = [-2 2];
            Y = [-2 2];
            Z = 0;
            conf.xref = [0 0 0];
            conf.secondary_sources.number = 56;
        case 'box'
            X = [-2 2];
            Y = [-2 2];
            Z = 0;
            conf.xref = [0 0 0];
            conf.secondary_sources.number = 80;
        case 'sphere'
            X = [-2 2];
            Y = 0;
            Z = [-2 2];
            conf.xref = [0 0 0];
            conf.secondary_sources.number = 900;
    end

    conf.secondary_sources.geometry = scenarios{ii,3};
    t = scenarios{ii,5};
    src = scenarios{ii,4};
    xs = scenarios{ii,5};
    % Adjust time for different source types (t=0 corresponds to first activity
    % of virtual source).
    if strcmp('ps',src) || strcmp('ls',src)
        t = 2/conf.c;   % time for traveling 2 m
    elseif strcmp('fs',src)
        t = 0.5/conf.c; % time for traveling 0.5 m
    else
        t = 0;
    end

    % ===== WFS ==========================================================
    if strcmp('WFS',scenarios{ii,1})
        % spatio-temporal impulse response
        try
            [p,x,y,z,x0] = sound_field_imp_wfs(X,Y,Z,xs,src,t,conf);
            if modus
                conf.plot.normalisation = 'max';
                plot_sound_field(p,X,Y,Z,x0,conf);
                title_str = sprintf('WFS %s %s array, %s, impulse response', ...
                    conf.dimension,conf.secondary_sources.geometry,src);
                title(title_str);
            end
        catch
            warning('%s: WFS impulse response %s array %s %s gives the following error message.', ...
                upper(mfilename),conf.secondary_sources.geometry,conf.dimension,src);
            lasterr
        end
    % ===== NFC-HOA ======================================================
    elseif strcmp('HOA',scenarios{ii,1})
        % spatio-temporal impulse response
%         try
            [p,x,y,z,x0] = sound_field_imp_nfchoa(X,Y,Z,xs,src,t,conf);
            if modus
                conf.plot.normalisation = 'max';
                plot_sound_field(p,X,Y,Z,x0,conf);
                title_str = sprintf('NFC-HOA %s %s array, %s, impulse response', ...
                    conf.dimension,conf.secondary_sources.geometry,src);
                title(title_str);
            end
%         catch
%             warning('%s: NFC-HOA impulse response %s array %s %s gives the following error message.', ...
%                 upper(mfilename),conf.secondary_sources.geometry,conf.dimension,src);
%             lasterr
%         end
    end
end


status = true;
