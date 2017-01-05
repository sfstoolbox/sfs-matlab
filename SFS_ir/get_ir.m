function [ir,x0] = get_ir(sofa,X,head_orientation,xs,coordinate_system,conf)
%GET_IR returns an impulse response for the given apparent angle
%
%   Usage: ir = get_ir(sofa,X,head_orientation,xs,[coordinate_system],conf)
%
%   Input parameters:
%       sofa              - impulse response data set (sofa struct/file)
%       X                 - position of the listener, specified in the defined
%                           coordinate_system, see below
%       head_orientation  - orientation of the listener with [phi theta] /
%                           (rad, rad)
%       xs                - position of the desired source, specified in the
%                           defined coordinate_system, see below. For
%                           SOFA convention SimpleFreeFieldHRIR xs will be
%                           interpreted relative to X, for MultiSpeakerBRIR as
%                           an absolute position.
%       coordinate_system - coordinate system X and xs are specified in,
%                           avialable systems are:
%                             'cartesian' - cartesian system (default) with
%                                           [x y z] / m
%                             'spherical' - spherical system with
%                                           [phi theta r] / (rad, rad, m)
%       conf              - configuration struct (see SFS_config)
%
%   Output parameters:
%       ir      - impulse response for the given position (length of IR x 2)
%       x0      - position corresponding to the returned impulse response
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
%             intpolate_ir, ir_correct_distance

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
    % desired distance. The desired direction is done by returning the nearest
    % neighbour or applying a linear interpolation.
    % NOTE: for SimpleFreeFieldHRIR head orientation is always zero in the SOFA
    % file and we handle a change in head orientation by changing the source
    % position accordingly.
    %
    % For SimpleFreeFieldHRIR only the relative position between listener
    % position and source position is of relevance.
    xs = xs-X+X_sofa;
    % Get measured loudspeaker positions and go to spherical coordinates
    x0 = sofa_get_secondary_sources(header,'spherical');
    [xs(1),xs(2),xs(3)] = cart2sph(xs(1),xs(2),xs(3));
    % Combine head orientation and desired direction of source (see note above)
    xs(1) = correct_azimuth(xs(1)-head_orientation(1));
    xs(2) = correct_elevation(xs(2)-head_orientation(2));
    % Find nearest neighbours and interpolate if desired and needed
    [neighbours,idx] = findnearestneighbour(x0(:,1:2)',xs(1:2),3);
    ir = sofa_get_data_fir(sofa,idx);
    ir = ir_correct_distance(ir,x0(idx,3),xs(3),conf);
    [ir,x0] = interpolate_ir(ir,neighbours,xs(1:2)',conf);
    [x0(1),x0(2),x0(3)] = sph2cart(x0(1),x0(2),xs(3));

elseif strcmp('MultiSpeakerBRIR',header.GLOBAL_SOFAConventions)
    %
    % http://www.sofaconventions.org/mediawiki/index.php/MultiSpeakerBRIR
    %
    % This looks for the nearest loudspeaker corresponding to the desired source
    % position. For that loudspeaker the impulse reponse for the specified head
    % orientation is returned. If the head orientation could not perfectly
    % matched, an interpolation is applied. If the specified head orientation is
    % out of bounds, the nearest head orientation is returned.
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
    [neighbours_emitter,idx_emitter] = findnearestneighbour(x0(:,1:3)',xs,1);
    x0 = neighbours_emitter(:,1);
    if norm(x0'-xs)>0.01
        warning('SFS:get_ir',['Your chosen loudspeaker position (%.2f,', ...
                              '%.2f,%.2f)m deviates from the measured ', ...
                              'one (%.2f,%.2f,%.2f)m.'], ...
                             xs(1),xs(2),xs(3),x0(1),x0(2),x0(3));
    end
    % Find nearest head orientation in the horizontal plane
    [phi,theta] = sofa_get_head_orientations(header);
    [x(:,1),x(:,2),x(:,3)] = sph2cart(phi,theta,ones(length(phi),1));
    [x_desired(1),x_desired(2),x_desired(3)] = ...
        sph2cart(head_orientation(1),head_orientation(2),1);
    [neighbours,idx_head] = findnearestneighbour(x',x_desired',3);
    [neighbours_head(1,:),neighbours_head(2,:)] = ...
        cart2sph(neighbours(1,:),neighbours(2,:),neighbours(3,:));
    % Check if head orientation is out of bounds
    if all(abs(head_orientation(1))>abs(neighbours_head(1,:)))
        warning('SFS:get_ir',['Head azimuth %.1f deg out of bound, ', ...
            'using %.1f deg instead.'], ...
            deg(head_orientation(1)), ...
            deg(neighbours_head(1,1)));
        head_orientation(1) = neighbours_head(1,1);
    end
    if all(abs(head_orientation(2))>abs(neighbours_head(2,:)))
        warning('SFS:get_ir',['Head elevation %.1f deg out of bound, ', ...
            'using %.1f deg instead.'], ...
            deg(head_orientation(2)), ...
            deg(neighbours_head(2,1)));
        head_orientation(2) = neighbours_head(2,1);
    end
    % Get the impulse responses, reshape and interpolate
    ir = sofa_get_data_fire(sofa,idx_head,idx_emitter);
    ir = reshape(ir,[size(ir,1) size(ir,2) size(ir,4)]); % [M R E N] => [M R N]
    ir = interpolate_ir(ir,neighbours_head,head_orientation',conf);

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
ir = squeeze(ir)';
% Convert x0 to the specified coordinate system
if strcmp('spherical',coordinate_system)
    [x0(1),x0(2),x0(3)] = cart2sph(x0(1),x0(2),x0(3));
end

warning('on','SOFA:upgrade')
