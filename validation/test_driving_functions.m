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
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************

% TODO: add mode to save data as reference data


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));


%% ===== Configuration ===================================================
conf = SFS_config;
L = 3;
f = 1000;
t = 200;
conf.useplot = modus;
conf.driving_functions = 'default';
% test scenarios
scenarios = { ...
    'WFS', '2D',   'linear',   'pw', [ 0.5  1.0  0.0]; ...
    'WFS', '2D',   'linear',   'ls', [ 0.0 -1.0  0.0]; ...
    'WFS', '2D',   'linear',   'ps', [ 0.0 -1.0  0.0]; ...
    'WFS', '2D',   'linear',   'fs', [ 0.0  1.0  0.0  0.0 1.0  0.0]; ...
    'WFS', '2D',   'circular', 'pw', [ 0.5  0.5  0.0]; ...
    'WFS', '2D',   'circular', 'ls', [ 0.0  2.5  0.0]; ...
    'WFS', '2D',   'circular', 'ps', [ 0.0  2.5  0.0]; ...
    'WFS', '2D',   'circular', 'fs', [ 0.0  0.5  0.0  0.0 -1.0  0.0]; ...
    'WFS', '2D',   'box',      'pw', [ 0.5  1.0  0.0]; ...
    'WFS', '2D',   'box',      'ls', [ 0.5  2.0  0.0]; ...
    'WFS', '2D',   'box',      'ps', [ 0.5  2.0  0.0]; ...
    'WFS', '2D',   'box',      'fs', [ 0.5  0.5  0.0  1.0  1.0  0.0]; ...
    'WFS', '2.5D', 'linear',   'pw', [ 0.5  1.0  0.0]; ...
    'WFS', '2.5D', 'linear',   'ps', [ 0.0 -1.0  0.0]; ...
    'WFS', '2.5D', 'linear',   'fs', [ 0.0  1.0  0.0  0.0 1.0  0.0]; ...
    'WFS', '2.5D', 'circular', 'pw', [ 0.5  0.5  0.0]; ...
    'WFS', '2.5D', 'circular', 'ps', [ 0.0  2.5  0.0]; ...
    'WFS', '2.5D', 'circular', 'fs', [ 0.0  0.5  0.0  0.0 -1.0  0.0]; ...
    'WFS', '2.5D', 'box',      'pw', [ 0.5  1.0  0.0]; ...
    'WFS', '2.5D', 'box',      'ps', [ 0.5  2.0  0.0]; ...
    'WFS', '2.5D', 'box',      'fs', [ 0.5  0.5  0.0  1.0  1.0  0.0]; ...
    %'HOA', '2.5D', 'circular', 'pw', [ 0.5  0.5  0.0]; ...
    %'HOA', '2.5D', 'circular', 'ps', [ 0.0  2.5  0.0]; ...
};

% Start testing
for ii=1:size(scenarios)

    % get current dimension
    conf.dimension = scenarios{ii,2};
   
    % get listening area
    switch scenarios{ii,3}
        case 'linear'
            X = [-2 2];
            Y = [-0.15 3];
            Z = 0;
            conf.xref = [0 2 0];
            conf.dx0 = 0.15;
        case 'circular'
            X = [-2 2];
            Y = [-2 2];
            Z = 0;
            conf.xref = [0 0 0];
            conf.dx0 = L*pi/56;
        case 'box'
            X = [-2 2];
            Y = [-2 2];
            Z = 0;
            conf.xref = [0 0 0];
            conf.dx0 = 0.15;
    end

    conf.array = scenarios{ii,3};
    src = scenarios{ii,4};
    xs = scenarios{ii,5};

    % ===== WFS ==========================================================
    if strcmp('WFS',scenarios{ii,1})
        % mono-frequent
        try
            [P,x,y,z,x0,win] = wave_field_mono_wfs(X,Y,Z,xs,src,f,L,conf);
            title_str = sprintf('WFS %s %s array, %s, mono-frequent', ...
                conf.dimension,conf.array,src);
            title(title_str);
        catch
            warning('%s: WFS mono-frequent %s array %s %s gives the following error message.', ...
                upper(mfilename),conf.array,conf.dimension,src);
            lasterr
        end
        % spatio-temporal impulse response
        try
            [p,x,y,z,x0,win] = wave_field_imp_wfs(X,Y,Z,xs,src,t,L,conf);
            title_str = sprintf('WFS %s %s array, %s, impulse response', ...
                conf.dimension,conf.array,src);
            title(title_str);
        catch
            warning('%s: WFS impulse response %s array %s %s gives the following error message.', ...
                upper(mfilename),conf.array,conf.dimension,src);
            lasterr
        end
    
    % ===== NFC-HOA ======================================================
    elseif strcmp('HOA',scenarios{ii,1})
        % mono-frequent
        try
            [P,x,y,z,x0,win] = wave_field_mono_nfchoa(X,Y,Z,xs,src,f,L,conf);
            title_str = sprintf('WFS %s %s array, %s, mono-frequent', ...
                conf.dimension,conf.array,src);
            title(title_str);
        catch
            warning('%s: WFS mono-frequent %s array %s %s gives the following error message.', ...
                upper(mfilename),conf.array,conf.dimension,src);
            lasterr
        end
        % spatio-temporal impulse response
        try
            [p,x,y,z,x0,win] = wave_field_imp_nfchoa(X,Y,Z,xs,src,t,L,conf);
            title_str = sprintf('WFS %s %s array, %s, impulse response', ...
                conf.dimension,conf.array,src);
            title(title_str);
        catch
            warning('%s: WFS impulse response %s array %s %s gives the following error message.', ...
                upper(mfilename),conf.array,conf.dimension,src);
            lasterr
        end
    end
end
