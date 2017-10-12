function status = test_interpolation_point_selection(modus)
%TEST_INTERPOLATION_POINT_SELECTION tests the correctness of the methods
%findconvexcone() and findvoronoi() for piecewise linear interpolation in
%3D grids
%
%   Usage: status = test_interpolation_point_selection(modus)
%
%   Input parameters:
%       modus   - 0: numerical (quiet)
%                 1: visual
%
%   Output parameters:
%       status - true or false

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2018 SFS Toolbox Developers                             *
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


%% ===== Secondary source grids ==========================================
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


%% ===== Definition of test cases ========================================
% NaN as reference is used as DON'T CARE and is not checked against the real
% value.
%
% --- Regular tests, should not return an error ---
regular_testcases{1} = {
    '3D grid, regular case: 3 points selected.', ...    % description
    'findconvexcone', ...                               % method
    x0_regular_3d, ...                                  % grid
    [2,0.1,0.1], ...                                    % desired point
    [9; 10; 1], ...                                     % reference indices
    [NaN; NaN; NaN]                                     % reference weights
    };
regular_testcases{2} = {
    '3D grid, regular case: 4 points selected, 1 has low weight.', ...
    'findvoronoi' ...
    x0_regular_3d, ...
    [2,0.1,0.1], ...
    [9; 10; 1; 32], ...
    [NaN; NaN; NaN; NaN]
    };

regular_testcases{3} = {
    '3D grid, xs coplanar with 2 points: 3 points selected, 1 has zero weight.', ...
    'findconvexcone', ...
    x0_regular_3d, ...
    [2,0.1,0], ...
    [9; 10; NaN], ... 
    [NaN; NaN; 0]
    };
regular_testcases{4} = {
    '3D grid, xs coplanar with 2 points: 4 points selected, 2 have low weight', ...
    'findvoronoi' ...
    x0_regular_3d, ...
    [2,0.1,0], ...
    [9; 10; 1; 33], ... 
    [NaN; NaN; NaN; NaN]
    };

regular_testcases{5} = {
    '3D grid, xs colinear with 1 point: 3 points selected, 2 have zero weight.', ...
    'findconvexcone', ...
    x0_regular_3d, ...
    [2,0,0], ...
    [9; NaN; NaN], ... 
    [NaN; 0; 0]
    };
regular_testcases{6} = {
    '3D grid, xs colinear with 1 point: 1 point selected.', ...
    'findvoronoi' ...
    x0_regular_3d, ...
    [2,0,0], ...
    [9], ... 
    [1]
    };

regular_testcases{7} = {
    '2D grid, regular case: 2 points selected.', ...
    'findconvexcone', ...
    x0_regular_2d, ...
    [2,0.1,0], ...
    [1; 2], ...
    [NaN; NaN]
    };
regular_testcases{8} = {
    '2D grid, regular case: 2 points selected.', ...
    'findvoronoi' ...
    x0_regular_2d, ...
    [2,0.1,0], ...
    [1; 2], ...
    [NaN; NaN]
    };

regular_testcases{9} = {
    '2D grid, xs colinear with 1 point: 2 points selected, 1 has zero weight.', ...
    'findconvexcone', ...
    x0_regular_2d, ...
    [2,0,0], ...
    [1; 24], ...
    [NaN; 0]
    };
regular_testcases{10} = {
    '2D grid, xs colinear with 1 point: 1 point selected.', ...
    'findvoronoi' ...
    x0_regular_2d, ...
    [2,0,0], ...
    [1], ...
    [1]
    };

regular_testcases{11} = {
    '2D grid, arc with gap smaller than 180 deg: 2 points selected.', ...
    'findconvexcone', ...
    x0_arc_2d, ...
    [2,-2,0], ...
    [1; 14], ...
    [NaN; NaN]
    };
regular_testcases{12} = {
    '2D grid, arc with gap smaller than 180 deg: 2 points selected.', ...
    'findvoronoi' ...
    x0_arc_2d, ...
    [2,-2,0], ...
    [1; 14], ...
    [NaN; NaN]
    };

regular_testcases{13} = {
    'Partial grid: 3 points selected, 1 has zero weight.', ...
    'findconvexcone', ...
    x0_upper_ring, ...
    [0,0,1], ...
    [NaN; NaN; NaN], ...
    [NaN; NaN; NaN]
    };
regular_testcases{14} = {
    'Partial grid: all points selected.', ...
    'findvoronoi' ...
    x0_upper_ring, ...
    [0,0,1], ...
    [NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN], ...
    [NaN; NaN; NaN; NaN; NaN; NaN; NaN; NaN]
    };

regular_testcases{15} = {
    'Partial grid, requested xs lies outside. Warning is issued.', ...
    'findconvexcone', ...
    x0_upper_ring, ...
    [1,0,0], ...
    [NaN; NaN], ...
    [NaN; NaN]
    };
regular_testcases{16} = {
    'Partial grid, requested xs lies outside. Warning is issued', ...
    'findvoronoi' ...
    x0_upper_ring, ...
    [1,0,0], ...
    [NaN; NaN; NaN], ...
    [NaN; NaN; NaN]
    };

regular_testcases{17} = {
    'Grid is not a sphere. Warning is issued.', ...
    'findconvexcone', ...
    x0_linear_shifted, ...
    [1,1,0], ...
    [NaN; NaN], ...
    [NaN; NaN]
    };
regular_testcases{18} = {
    'Grid is not a sphere. Warning is issued', ...
    'findvoronoi' ...
    x0_linear_shifted, ...
    [1,1,0], ...
    [10], ...
    [1]
    };

% --- Erroneous tests, should return an error ---
erroneous_testcases{1} = {
    'Grid and xs are colinear. convhulln raises error', ...
    'findconvexcone', ...
    x0_linear, ...
    [0.1,0,0]
    };
erroneous_testcases{2} = {
    'Grid and xs are colinear. convhulln raises error', ...
    'findvoronoi', ...
    x0_linear, ...
    [0.1,0,0]
    };


%% ===== Run tests =======================================================
for testcase_tmp = regular_testcases
    testcase = testcase_tmp{1};
    desc_str = testcase{1};
    method = testcase{2};
    x0 = testcase{3};
    xs = testcase{4};
    x0_indices_ref = testcase{5};
    x0_weights_ref = testcase{6};
    if modus
        disp(['method: ', method,  ', test case: ', desc_str]);
    end
    
    if strcmp('findvoronoi',method)
        [indices,weights] = findvoronoi(x0,xs);
    elseif strcmp('findconvexcone',method)
        [indices,weights] = findconvexcone(x0,xs);
    end
   
    if modus
        plot_point_selection(x0,xs,indices,weights,desc_str,method);
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

for testcase_tmp = erroneous_testcases
    testcase = testcase_tmp{1};
    desc_str = testcase{1};
    method = testcase{2};
    x0 = testcase{3};
    xs = testcase{4};
    if modus
        plot_point_selection(x0,xs,[],[],desc_str,method);
        disp(['method: ', method, ', test case: ' , desc_str]);
    end
    
    if strcmp('findvoronoi',method)
        try
            findvoronoi(x0,xs)
        catch
            continue
        end
    elseif strcmp('findconvexcone',method)
        try
            findconvexcone(x0,xs)
        catch
            continue
        end
    end
end


status = true;


end


%% ===== Functions =======================================================
function plot_point_selection(x0,xs,indices,weights,desc_str,method)
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
    legend('x0','xs','weights','selected points','center' );
    title({desc_str; ['Method: ', method]},'interpreter','none');
end
