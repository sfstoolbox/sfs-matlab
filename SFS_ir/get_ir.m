function ir = get_ir(sofa,X,head_orientation,xs,coordinate_system,conf)
%GET_IR returns an impulse response for the given apparent angle
%
%   Usage: ir = get_ir(sofa,X,head_orientation,xs,coordinate_system,conf)
%
%   Input parameters:
%       sofa               - impulse response data set (sofa struct/file)
%       X                  - position of the listener, specified in the defined
%                            coordinate_system, see below
%       head_orientation   - orientation of the listener with [phi theta] /
%                            (rad, rad)
%       xs                 - position of the desired source, specified in the
%                            defined coordinate_system, see below. For
%                            SOFA convention SimpleFreeFieldHRIR xs will be
%                            interpreted relative to X, for MultiSpeakerBRIR as
%                            an absolute position.
%       coordinate_system  - coordinate system X and xs are specified in,
%                            avialable systems are:
%                              'cartesian' - cartesian system (default) with
%                                            [x y z] / m
%                              'spherical' - spherical system with
%                                            [phi theta r] / (rad, rad, m)
%       conf               - configuration struct (see SFS_config)
%
%   Output parameters:
%       ir      - impulse response for the given position (length of IR x 2)
%
%   GET_IR(sofa,X,head_orientation,xs,conf) returns a single impulse response
%   from the given SOFA file or struct. The impulse response is determined by
%   the position X and head orientation head_orientation of the listener, and
%   the position xs of the desired point source.
%   If the desired angles are not present in the SOFA data set and
%   conf.ir.useinterpolation is set to true an interpolation is applied to
%   create the impulse response. For the SOFA convention MultiSpeakerBRIR the
%   interpolation is only performed for the head orientations not the different
%   loudspeaker positions.
%
%   For the SOFA convention SimpleFreeFieldHRIR the desired distance between
%   the point source and listener is achieved by delaying and weighting the
%   impulse response. In this case zeros are padded at the beginning and end
%   of the impulse response and the length of ir is given by conf.N. The zeros
%   added at the beginning correspond to the actual radius of the measured
%   impulse response.
%
%   For a description of the SOFA file format see: http://sofaconventions.org
%
%   See also: SOFAload, sofa_get_header, sofa_get_data_fir, sofa_get_data_fire,
%             interpolate_ir, ir_correct_distance

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
nargmin = 5;
nargmax = 6;
narginchk(nargmin,nargmax)
if nargin<nargmax
    conf = coordinate_system;
    coordinate_system = 'cartesian';
end
if length(head_orientation)==1
    head_orientation = [head_orientation 0]; % [azimuth elevation]
end


%% ===== Computation ====================================================
warning('off','SOFA:upgrade')
%
% === SOFA meta data ===
header = sofa_get_header(sofa);

% === Coordinate system conversion ===
% Convert everything to cartesian coordinates
if strcmp('spherical',coordinate_system)
    [xs(1),xs(2),xs(3)] = sph2cart(xs(1),xs(2),xs(3));
    [X(1),X(2),X(3)] = sph2cart(X(1),X(2),X(3));
elseif ~strcmp('cartesian',coordinate_system)
    error('%s: unknown coordinate system type.',upper(mfilename));
end
% Get the listener position during the measurement
X_sofa = sofa_get_listener_position(header,'cartesian');

% === Get Impulse Response ===
if strcmp('SimpleFreeFieldHRIR',header.GLOBAL_SOFAConventions)
    %
    % http://www.sofaconventions.org/mediawiki/index.php/SimpleFreeFieldHRIR
    %
    % Returns a single impulse response for the desired position. The impulse
    % response is shifted in time and its amplitude is weighted according to the
    % desired distance. The desired direction is handled by returning the nearest
    % neighbour or applying a linear interpolation.
    % NOTE: for SimpleFreeFieldHRIR head orientation is always zero in the SOFA
    % file and we handle a change in head orientation by changing the source
    % position accordingly.
    %
    % For SimpleFreeFieldHRIR only the relative position between the listener
    % and the source position is of relevance.
    xs = xs-X;
    % Include head orientation of listener into relative position
    [xs(1),xs(2),xs(3)] = cart2sph(xs(1),xs(2),xs(3));
    xs(1) = correct_azimuth(xs(1)-head_orientation(1));
    xs(2) = correct_elevation(xs(2)-head_orientation(2));
    [xs(1),xs(2),xs(3)] = sph2cart(xs(1),xs(2),xs(3));
    % Get measured loudspeaker positions relative to dummy head position
    x0 = sofa_get_secondary_sources(header,'cartesian');
    x0 = bsxfun(@minus,x0(:,1:3),X_sofa);
    % Get nearest neighbour point or points for interpolation
    [idx,weights] = point_selection(x0,xs,conf);
    % Get impulse responses according to point selection
    ir = sofa_get_data_fir(sofa,idx);
    % Correct distance
    [x0(:,1),x0(:,2),x0(:,3)] = cart2sph(x0(:,1),x0(:,2),x0(:,3));
    [xs(1),xs(2),xs(3)] = cart2sph(xs(1),xs(2),xs(3));
    ir = ir_correct_distance(ir,x0(idx,3),xs(3),conf);
    [x0(:,1),x0(:,2),x0(:,3)] = sph2cart(x0(:,1),x0(:,2),x0(:,3));
    % Select or interpolate to desired impulse response
    ir = interpolate_ir(ir,weights,conf);

elseif strcmp('MultiSpeakerBRIR',header.GLOBAL_SOFAConventions)
    %
    % http://www.sofaconventions.org/mediawiki/index.php/MultiSpeakerBRIR
    %
    % This looks for the nearest loudspeaker corresponding to the desired source
    % position. For that loudspeaker, the impulse reponse for the specified head
    % orientation is returned. If the head orientation could not be perfectly
    % matched, an interpolation is applied if desired.
    %
    % Check if we are requesting a listening position that is available
    if norm(X-X_sofa)>0.01
        warning('SFS:get_ir',['Your chosen listening position (%.2f,', ...
                              '%.2f,%.2f)m is not available, using the ', ...
                              'measured one (%.2f,%.2f,%.2f)m instead.'], ...
                             X(1),X(2),X(3),X_sofa(1),X_sofa(2),X_sofa(3));
    end
    % Find nearest loudspeaker
    x0 = sofa_get_secondary_sources(header,'cartesian');
    x0 = x0(:,1:3); % need position only
    idx_emitter = findnearestneighbour(x0,xs,1);
    x0_new = x0(idx_emitter,:);
    if norm(x0_new-xs)>0.01
        warning('SFS:get_ir',['Your chosen loudspeaker position (%.2f,', ...
                              '%.2f,%.2f)m deviates from the measured ', ...
                              'one (%.2f,%.2f,%.2f)m.'], ...
                             xs(1),xs(2),xs(3),x0_new(1),x0_new(2),x0_new(3));
    end
    % Get nearest neighbour point or points for interpolation for head orientation
    [phi,theta,r] = sofa_get_head_orientations(header);
    [sofa_head_orientations(:,1),sofa_head_orientations(:,2),...
        sofa_head_orientations(:,3)] = sph2cart(phi,theta,r);
    [head_orientation(1),head_orientation(2),head_orientation(3)] = ...
        sph2cart(head_orientation(1),head_orientation(2),1);
    [idx_head,weights] = point_selection(...
        sofa_head_orientations,head_orientation,conf);
    % Get impulse responses according to point selection
    ir = sofa_get_data_fire(sofa,idx_head,idx_emitter);
    ir = reshape(ir,[size(ir,1) size(ir,2) size(ir,4)]); % [M R E N] => [M R N]
    % Select or interpolate to desired impulse response
    ir = interpolate_ir(ir,weights,conf);

elseif strcmp('SingleRoomDRIR',header.GLOBAL_SOFAConventions)
    %
    % http://www.sofaconventions.org/mediawiki/index.php/SingleRoomDRIR
    %
    error(['%s: SingleRoomDRIR is not supported as it should handle ', ...
           'microphone array recordings. If you used it for (multiple) ', ...
           'loudspeakers in a room you should consider to use ', ...
           'MultiSpeakerBRIR instead.'], upper(mfilename));

else
    error('%s: %s convention is currently not supported.', ...
        upper(mfilename),header.GLOBAL_SOFAConventions);
end

% Reshape [1 2 N] to [N 2]
ir = squeeze(ir).';

warning('on','SOFA:upgrade')
end

% =========================================================================

function [idx,weights] = point_selection(x0,xs,conf)
%POINT_SELECTION selects points out of x0 with weights that can be used to
% interpolate an impulse response associated with xs out of impulse responses
% associated with x0(idx,:). Instead of interpolation, the nearest neighbour point
% can be selected as well.
%
%   Input parameters:
%       x0       - points in R^3 / m [nx3]
%       xs       - desired point in R^3 / m [1x3]
%       conf     - 
%
%   Output parameters:
%       idx      - row indices of N points in A [Nx1]
%       weights  - weights [Nx1]
useinterpolation = conf.ir.useinterpolation;
% Check for old configuration
if useinterpolation
    if ~isfield(conf.ir,'interpolationpointselection')
        warning('SFS:get_ir:interpolationpointselection',...
            ['%s: no method for selection of interpolation points provided, ', ...
            'will use method ''nearestneighbour''.'],upper(mfilename));
        interpolationpointselection = 'nearestneighbour';
    else
        interpolationpointselection = conf.ir.interpolationpointselection;
    end
end
% Get nearest neighbour point or points for interpolation
if useinterpolation
    switch interpolationpointselection
    case 'nearestneighbour'
        [idx,weights] = findnearestneighbour(x0,xs,2);
    case 'delaunay'
        [idx,weights] = findconvexcone(x0,xs);
    case 'voronoi'
        to_be_implemented;
    otherwise
        error(['%s: ''%s'' is an unknown method to select interpolation ', ...
            'points.'],upper(mfilename),interpolationpointselection);
    end
else
    [idx,weights] = findnearestneighbour(x0,xs,1);
end
end
