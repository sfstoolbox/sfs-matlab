function status = test_secondary_source_selection(modus)
%TEST_SECONDARY_SOURCE_SELECTION tests the correctness of
%secondary_source_selection()
%
%   Usage: status = test_secondary_source_selection(modus)
%
%   Input parameters:
%       modus   - 0: numerical (quiet)
%                 1: visual (not available)
%                 2: numerical verbose
%
%   Output parameters:
%       status - true or false
%
%   TEST_SECONDARY_SOURCE_SELECTION(modus) checks if the secondary source
%   selection needed for Wave Field Syntesis is implemented correctly in
%   secondary_source_selection().

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
        return;
    end
elseif modus==2
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
end


status = true;
