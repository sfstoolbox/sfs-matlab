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

% Create some grids
[X,Y,Z] = sph2cart((0:pi/4:7/4*pi).',pi/4*ones(8,1),ones(8,1));
x0_upper_ring = [X,Y,Z];
[X,Y,Z] = sph2cart((0:pi/12:23/12*pi).',zeros(24,1),ones(24,1));
x0_center_ring = [X,Y,Z];

x0_linear = [linspace(-1,1,10).',zeros(10,2)];
x0_linear_shifted = [linspace(-1,1,10).',ones(10,1),zeros(10,1)];
south_pole = [0,0,-1];

x0_regular_3d = [x0_upper_ring; x0_center_ring; south_pole];
x0_regular_2d = x0_center_ring*rotation_matrix(pi/100,1);

x0_arc_2d = x0_regular_2d(1:14,:);


% Test cases:
% NaN as reference is used as DON'T CARE
% (NaN as result shall not occur)
regular_testcases{1}= {
    '3D grid: regular case: 3 points selected.', ... % 1. description
    x0_regular_3d, ... % 2. grid
    [2,0.1,0.1], ...  % 3. desired point
    [9; 10; 1], ...   % 4. reference indices
    [NaN; NaN; NaN]   % 5. reference weights
    };
regular_testcases{2} = {
    '3D grid: xs coplanar with 2 points: 3 points selected, 1 has zero weight.', ...
    x0_regular_3d, ...
    [2,0.1,0], ...
    [9; 10; NaN], ... 
    [NaN; NaN; 0]
    };
regular_testcases{3}= {
    '3D grid: xs colinear with 1 point:  3 points selected, 2 have zero weight.', ...
    x0_regular_3d, ...
    [2,0,0], ...
    [9; NaN; NaN], ... 
    [NaN; 0; 0]
    };
regular_testcases{4}= {
    '2D grid: regular case: 2 points selected.', ...
    x0_regular_2d, ...
    [2,0.1,0], ...
    [1; 2], ...
    [NaN; NaN]
    };
regular_testcases{5}= {
    '2D grid: xs colinear with 1 point:  2 points selected, 1 has zero weight.', ...
    x0_regular_2d, ...
    [2,0,0], ...
    [1; 24], ...
    [NaN; 0]
    };
regular_testcases{6}= {
    '2D grid: arc with gap smaller than 180 deg: 2 points selected.', ...
    x0_arc_2d, ...
    [2,-2,0], ...
    [1; 14], ...
    [NaN; NaN]
    };
regular_testcases{7}= {
    'Partial grid.', ...
    x0_upper_ring, ...
    [0,0,1], ...
    [NaN; NaN; NaN], ...
    [NaN; NaN; NaN]
    };
regular_testcases{8}= {
    'Partial grid, requested xs lies outside. Warning is issued.', ...
    x0_upper_ring, ...
    [1,0,0], ...
    [NaN; NaN], ...
    [NaN; NaN]
    };
regular_testcases{9}= {
    'Grid is not a sphere. Warning is issued.', ...
    x0_linear_shifted, ...
    [1,1,0], ...
    [NaN; NaN], ...
    [NaN; NaN]
    };

for testcase_tmp = regular_testcases
    testcase = testcase_tmp{1};
    desc_str= testcase{1};
    x0 = testcase{2};
    xs = testcase{3};
    x0_indices_ref = testcase{4};
    x0_weights_ref = testcase{5};
    if modus
        disp(['test case: ' , desc_str]);
    end
    [indices,weights] = findconvexcone(x0,xs);
    if modus
        plot_point_selection(x0,xs,indices,weights,desc_str);
    end
    if any(weights < 0)
    error('%s: In %s: negative weights. ', ...
        upper(mfilename),desc_str);
    end
    if ~(all(indices == x0_indices_ref | isnan(x0_indices_ref)))
    error('%s: In %s: wrong indices ', ...
        upper(mfilename),desc_str);
    end
    if ~(all(weights == x0_weights_ref| isnan(x0_weights_ref)))
    error('%s: In %s: wrong weights ', ...
        upper(mfilename),desc_str);
    end
end


erroneous_testcases{1}= {
    'x0 colinear through origin. convhulln raises error', ...
    x0_linear, ...
    [0.1,0,0]
    };

for testcase_tmp = erroneous_testcases
    testcase = testcase_tmp{1};
    desc_str= testcase{1};
    x0 = testcase{2};
    xs = testcase{3};
    if modus
        plot_point_selection(x0,xs,[],[],desc_str);
        disp(['test case: ' , desc_str]);
    end
    try
        findconvexcone(x0,xs)
        return
    end
end

status = true;

    function plot_point_selection(x0,xs,indices,weights,desc_str)
        if isoctave
            point_size = 12;
            selected_point_size = weights*20 + 1;
        else
            point_size = 100;
            selected_point_size = weights*100 + 1;
        end
        figure
        scatter3(x0(:,1),x0(:,2),x0(:,3),point_size,'b','.');
        hold on
        quiver3(0,0,0,xs(1),xs(2),xs(3),'k');
        scatter3(x0(indices,1),x0(indices,2),x0(indices,3),selected_point_size,'r');
        scatter3(x0(indices,1),x0(indices,2),x0(indices,3),point_size,'r','x');
        scatter3(0,0,0,point_size,'k','.');
        hold off
        axis equal
        xlabel('x');
        ylabel('y');
        zlabel('z');
        title(desc_str,'interpreter','none');
    end
end
