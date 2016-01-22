function boolean = test_secondary_source_selection(modus)
%TEST_SECONDARY_SOURCE_SELECTION tests the correctness of
%secondary_source_selection()
%
%   Usage: boolean = test_secondary_source_selection(modus)
%
%   Input parameters:
%       modus   - 0: numerical (quiet)
%                 1: numerical (verbose)
%                 2: visual
%
%   Output parameters:
%       booelan - true or false
%
%   TEST_SECONDARY_SELECTION(modus) checks if the secondary source selection
%   needed for Wave Field Syntesis is implemented correctly in
%   secondary_source_selection().

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
ref_selection_circular_pw = ...
    [0   1   1   1   1   1   1   1   0   0   0   0   0   0   0   0]';
ref_selection_circular_ps = ...
    [0   0   1   1   1   1   1   0   0   0   0   0   0   0   0   0]';
ref_selection_circular_fs = ...
    [1   1   1   1   1   0   0   0   0   0   0   0   0   0   0   0]';
ref_selection_linear_pw = ...
    [1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1]';
ref_selection_linear_ps = ...
    [1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1]';
ref_selection_linear_fs = ...
    [0   0   0   0   0   0   0   0   1   1   1   1   1   1   1   1]';
ref_selection_box_pw = ...
    [0   0   0   0   1   1   1   1   0   0   0   0   0   0   0   0]';
ref_selection_box_ps = ...
    [0   0   0   0   1   1   1   1   0   0   0   0   0   0   0   0]';
ref_selection_box_fs = ...
    [0   1   1   1   1   1   1   0   0   0   0   0   0   0   0   0]';
% Calculate current values
% circular array
conf.secondary_source.geometry = 'circular';
x0 = secondary_source_positions(conf);
[~,selection_circular_pw] = secondary_source_selection(x0,[0 -1 0],'pw');
[~,selection_circular_ps] = secondary_source_selection(x0,[0 2.5 0],'ps');
[~,selection_circular_fs] = secondary_source_selection(x0,[0.5 0.5 0 -1 -1 0],'fs');
% linear array
conf.secondary_sources.geometry = 'linear';
x0 = secondary_source_positions(conf);
[~,selection_linear_pw] = secondary_source_selection(x0,[0 -1 0],'pw');
[~,selection_linear_ps] = secondary_source_selection(x0,[0 1 0],'ps');
[~,selection_linear_fs] = secondary_source_selection(x0,[0.5 -0.5 0 -1 -1 0],'fs');
% box form array
conf.secondary_sources.geometry = 'box';
x0 = secondary_source_positions(conf);
[~,selection_box_pw] = secondary_source_selection(x0,[0 -1 0],'pw');
[~,selection_box_ps] = secondary_source_selection(x0,[0 3.5 0],'ps');
[~,selection_box_fs] = secondary_source_selection(x0,[0.5 0.5 0 -1 -1 0],'fs');

if modus==0
    % Numerical mode (quiet)
    if ~all(eq(ref_selection_circular_pw,selection_circular_pw)) || ...
       ~all(eq(ref_selection_circular_ps,selection_circular_ps)) || ...
       ~all(eq(ref_selection_circular_fs,selection_circular_fs)) || ...
       ~all(eq(ref_selection_linear_pw,selection_linear_pw)) || ...
       ~all(eq(ref_selection_linear_ps,selection_linear_ps)) || ...
       ~all(eq(ref_selection_linear_fs,selection_linear_fs)) || ...
       ~all(eq(ref_selection_box_pw,selection_box_pw)) || ...
       ~all(eq(ref_selection_box_ps,selection_box_ps)) || ...
       ~all(eq(ref_selection_box_fs,selection_box_fs))
        boolean = false;
    end
elseif modus==1
    message = 'wrong secondary source selection for a';
    if ~all(eq(ref_selection_circular_pw,selection_circular_pw))
        error('%s: %s circular array and a plane wave.', ...
            upper(mfilename),message);
    end
    if ~all(eq(ref_selection_circular_ps,selection_circular_ps))
        error('%s: %s circular array and a point source.', ...
            upper(mfilename),message);
    end
    if ~all(eq(ref_selection_circular_fs,selection_circular_fs))
        error('%s: %s circular array and a focused source.', ...
            upper(mfilename),message);
    end
    if ~all(eq(ref_selection_linear_pw,selection_linear_pw))
        error('%s: %s linear array and a plane wave.', ...
            upper(mfilename),message);
    end
    if ~all(eq(ref_selection_linear_ps,selection_linear_ps))
        error('%s: %s linear array and a point source.', ...
            upper(mfilename),message);
    end
    if ~all(eq(ref_selection_linear_fs,selection_linear_fs))
        error('%s: %s linear array anda focused source.', ...
            upper(mfilename),message);
    end
    if ~all(eq(ref_selection_box_pw,selection_box_pw))
        error('%s: %s box shaped array and a plane wave.', ...
            upper(mfilename),message);
    end
    if ~all(eq(ref_selection_box_ps,selection_box_ps))
        error('%s: %s box shaped array and a point source.', ...
            upper(mfilename),message);
    end
    if ~all(eq(ref_selection_box_fs,selection_box_fs))
        error('%s: %s box shaped array and a focused source.', ...
            upper(mfilename),message);
    end
else
    error('%s: modus has to be 0 (numerical quiet), 1 (numerical), ', ...
          upper(mfilename));
end
