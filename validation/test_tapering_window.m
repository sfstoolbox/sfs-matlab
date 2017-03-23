function status = test_tapering_window(modus)
%TEST_TAPERING_WINDOW tests the tapering_window() function for applying tapering
%to the secondary sources in Wave Field Synthesis
%
%   Usage: status = test_tapering_window(modus)
%
%   Input parameters:
%       modus   - 0: numerical (quiet)
%                 1: visual (not available)
%                 2: numerical verbose
%
%   Output parameters:
%       status - true or false
%
%   TEST_TAPERING_WINDOW(modus) checks if the tapering window applied to the
%   secondary sources in Wave Field Synthesis is working.

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


%% ===== Main ============================================================
conf = SFS_config;
conf.secondary_sources.number = 16;
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
        return;
    end
elseif modus==2
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


status = true;
