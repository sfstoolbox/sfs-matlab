function [diam,center] = secondary_source_diameter(conf)
%SECONDARY_SOURCE_DIAMETER calculates the maximum distance
% between the secondary sources (the diameter) and the center of the
% smallest ball that contains the array.
%
%   Usage: [diam,center] = secondary_source_diameter(conf)
%
%   Input parameters:
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       diam        - diameter of secondary source distribution / m
%       center      - center of the ball containing SSD / m [1x3]
%
%   SECONDARAY_SOURCE_DIAMETER(conf) calculates the maximum
%   Euklidian distance between the given secondary sources. Additionaly,
%   the center of the encompassing ball is returned. If one of the predefined
%   secondary source distributions 'linear', 'circular', or 'spherical' is used,
%   the returned diameter is equal to conf.secondary_sources.size.
%   If 'box' or 'rounded-box' is used, the diameter of the bounding box
%   is returned.
%
%   See also: driving_function_imp_wfs, secondary_source_positions

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


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargstruct(conf);


%% ===== Configuration ==================================================
geometry = conf.secondary_sources.geometry;


%% ===== Calculation ====================================================
if strcmp('line',geometry)   || strcmp('linear',geometry)    || ...
   strcmp('circle',geometry) || strcmp('circular',geometry)  || ...
   strcmp('sphere',geometry) || strcmp('spherical',geometry)
    diam = conf.secondary_sources.size;
    center = conf.secondary_sources.center;
elseif strcmp('box',geometry)
    dx0 = conf.secondary_sources.size/(conf.secondary_sources.number/4-1);
    diam = (conf.secondary_sources.size+dx0)*sqrt(2);
    center = conf.secondary_sources.center;
elseif strcmp('rounded-box',geometry)
    diam = (conf.secondary_sources.size)*sqrt(2);
    center = conf.secondary_sources.center;
else
    x0 = conf.secondary_sources.x0;
    if isempty(x0)
        error(['%s: conf.secondary_sources.x0 must contain the secondary ',...
            'sources when using geometry %s.'],upper(mfilename),geometry);
    end
    % Find source1 :=  source with largest distance from origin
    [~,idx1] = max(vector_norm(x0(:,1:3),2));
    % Find source2 := source with maximum distance to source1
    [diam,idx2] = max(vector_norm(x0(:,1:3) - ...
        repmat(x0(idx1,1:3),[size(x0,1),1]),2));
    % Center is half-way between source1 and source2
    center = x0(idx1,1:3) +  0.5 * (x0(idx2,1:3) - x0(idx1,1:3));
end
