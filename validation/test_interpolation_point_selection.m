function status = test_interpolation_point_selection(modus)
%TEST_INTERPOLATION_POINT_SELECTION tests the correctness of
%findconvexcone() for piecewise linear interpolation in 3D grids
%
%   Usage: status = test_interpolation_point_selection(modus)
%
%   Input parameters:
%       modus   - 0: numerical (quiet)
%                 1: visual
%
%   Output parameters:
%       status - true or false
%
%   TEST_SECONDARY_SOURCE_SELECTION(modus) checks if the grid point
%   selection in findconvexcone for piecewise linear interpolation 
%   is implemented correctly

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Developers                             *
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

% create grid: two rings with differenent elevation
[X,Y,Z] = sph2cart((0:pi/4:7/4*pi).',pi/4*ones(8,1),ones(8,1));
x0_two_rings = [X,Y,Z];
[X,Y,Z] = sph2cart((0:pi/12:23/12*pi).',zeros(24,1),ones(24,1));
x0_two_rings = [x0_two_rings; [X,Y,Z]];

% Test cases:
% NaN as reference is used as DON'T CARE
% (NaN as result shall not occur)
testcases{1}= {
    'regular case: 3 points selected', ... % 1. description
    x0_two_rings, ... % 2. grid
    [2,0.1,0.1], ...  % 3. desired point
    [9; 10; 1], ...   % 4. reference indeces
    [NaN; NaN; NaN]   % 5. reference weights
    };
testcases{2} = {
    'degenerate case: xs colinear with 2 points', ...
    x0_two_rings, ...
    [2,0.1,0], ...
    [9; 10; NaN], ... 
    [NaN; NaN; 0]
    };
testcases{3}= {
    'degenerate case: xs coincident', ...
    x0_two_rings, ...
    [2,0,0], ...
    [9; NaN; NaN], ... 
    [NaN; 0; 0]
    };


for testcase_tmp = testcases
    testcase = testcase_tmp{1};
    desc_str= testcase{1};
    x0 = testcase{2};
    xs = testcase{3};
    x0_indeces_ref = testcase{4};
    x0_weights_ref = testcase{5};
    [indeces,weights] = findconvexcone(x0,xs);
    if modus
        plot_point_selection(x0,xs,indeces,weights,desc_str);
    end
    if any(weights < 0)
    error('%s: In %s: negative weights. ', ...
        upper(mfilename),desc_str);
    end
    if ~(all(indeces == x0_indeces_ref | isnan(x0_indeces_ref)))
    error('%s: In %s: wrong indeces ', ...
        upper(mfilename),desc_str);
    end
    if ~(all(weights == x0_weights_ref| isnan(x0_weights_ref)))
    error('%s: In %s: wrong weights ', ...
        upper(mfilename),desc_str);
    end
end

status = true;
    
    function plot_point_selection(x0,xs,indeces,weights,desc_str)
        point_size = weights*100 + 1;
        figure
        scatter3(x0(:,1),x0(:,2),x0(:,3),'b.');
        hold on
        quiver3(0,0,0,xs(1),xs(2),xs(3),'k');
        scatter3(x0(indeces,1),x0(indeces,2),x0(indeces,3),point_size);
        hold off
        axis equal
        xlabel('x');
        ylabel('y');
        zlabel('z');
        title(desc_str);
    end
end
