function boolean = test_tapering_window(modus)
%TEST_TAPERING_WINDOW tests the tapering_window() function for applying tapering
%to the secondary sources in Wave Field Synthesis
%
%   Usage: boolean = test_tapering_window(modus)
%
%   Input parameters:
%       modus   - 0: numerical (quiet)
%                 1: numerical (verbose)
%
%   Output parameters:
%       booelan - true or false
%
%   TEST_TAPERING_WINDOW(modus) checks if the tapering window applied to the
%   secondary sources in Wave Field Synthesis is working.

%*****************************************************************************
% Copyright (c) 2010-2016 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2016 Institut fuer Nachrichtentechnik                   *
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
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);


%% ===== Main ============================================================
conf = SFS_config;
conf.secondary_sources.number = 16;
boolean = true;
% reference values
ref_win_linear = [ 
   0.25000
   1.00000
   1.00000
   1.00000
   1.00000
   1.00000
   1.00000
   1.00000
   1.00000
   1.00000
   1.00000
   1.00000
   1.00000
   1.00000
   1.00000
   0.25000
];
ref_win_box = [
   1
   1
   1
   1
   1
   1
   1
   1
];
ref_win_circular1 = [
   1
   1
   1
   1
   1
];
ref_win_circular2 = [
   1
   1
   1
   1
   1
];
% Calculate current values
% linear secondary sources
conf.secondary_sources.geometry = 'linear';
x0 = secondary_source_positions(conf);
win_linear = tapering_window(x0,conf);
% box shaped secondary sources
conf.secondary_sources.geometry = 'box';
x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,[-1 -1 0],'pw');
win_box = tapering_window(x0,conf);
% circular secondary sources
conf.secondary_sources.geometry = 'circular';
x0 = secondary_source_positions(conf);
x01 = secondary_source_selection(x0,[0 2.5 0],'ps');
win_circular1 = tapering_window(x01,conf);
x02 = secondary_source_selection(x0,[2.5 0 0],'ps');
win_circular2 = tapering_window(x02,conf);

if modus==0
    % Numerical mode (quiet)
    if sum(abs(ref_win_circular1-win_circular1))>eps || ...
       sum(abs(ref_win_circular2-win_circular2))>eps || ...
       sum(abs(ref_win_linear-win_linear))>eps || ...
       sum(abs(ref_win_box-win_box))>eps
        boolean = false;
    end
elseif modus==1
    message = 'wrong tapering window for';
    if sum(abs(ref_win_circular1-win_circular1))>eps || ...
       sum(abs(ref_win_circular2-win_circular2))>eps
        error('%s: %s circular secondary sources.', ...
            upper(mfilename),message);
    end
    if sum(abs(ref_win_linear-win_linear))>eps
        error('%s: %s linear secondary sources.', ...
            upper(mfilename),message);
    end
    if sum(abs(ref_win_box-win_box))>eps
        error('%s: %s box shaped secondary sources.', ...
            upper(mfilename),message);
    end
else
    error('%s: modus has to be 0 (numerical quiet), 1 (numerical), ', ...
          upper(mfilename));
end
