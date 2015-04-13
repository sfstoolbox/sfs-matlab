function ir = get_ir(sofa,X,head_orientation,xs,coordinate_system,conf)
%GET_IR returns an impulse response for the given apparent angle
%
%   Usage: ir = get_ir(sofa,X,head_orientation,xs,[coordinate_system],[conf])
%
%   Input parameters:
%       sofa              - impulse response data set (sofa struct/file)
%       X                 - position of the listener, specified in the defined
%                           coordinate_system, see below
%       head_orientation  - orientation of the listener with [phi theta] /
%                           (rad, rad)
%       xs                - position of the desired source, specified in the
%                           defined coordinate_system, see below
%       coordinate_system - coordinate system xs is specified in, avialable
%                           systems are:
%                             'spherical' - spherical system (default) with
%                                           [phi theta r] / (rad, rad, m)
%                             'cartesian' - cartesian system with [x y z] / m
%       conf              - optional configuration struct (see SFS_config),
%                           which will be passed to:
%                             interpolate_ir
%                             ir_correct_distance (only for SimpleFreeFieldHRIR)
%
%   Output parameters:
%       ir      - impulse response for the given position (length of IR x 2)
%
%   GET_IR(sofa,X,head_orientation,xs) returns a single impulse response from
%   the given SOFA file or struct. The impulse response is determined by the
%   position X and head orientation head_orientation of the listener, and the
%   position xs of the desired point source.
%   For the SOFA convention SimpleFreeFieldHRIR the desired distance between
%   the point source and listener is achieved by delaying and weighting the
%   impulse response. Distances larger than 10m are ignored and set constantly
%   to 10m.
%   If the desired angles are not present in the SOFA data set and
%   conf.ir.useinterpolation is set to true an interpolation is applied to
%   create the impulse response. For the SOFA convention MultiSpeakerBRIR the
%   interpolation is only performed for the head orientations not the different
%   loudspeaker positions.
%   A further configuration setting that is considered is conf.ir.useoriglength,
%   which indicates if additional zeros corresponding to the actual radius of
%   the measured impulse responses should be added at the beginning of every
%   impulse response (if set to false). If you know that the measured impulse
%   responses include already the zeros from the measurement it can be safely
%   set to true. This is important because the delaying of the impulse responses
%   in order to achieve the correct distance require enough zeros at the
%   beginning of every impulse response.
%
%   For a description of the SOFA file format see: http://sofaconventions.org
%
%   See also: SOFAload, sofa_get_header, sofa_get_data_fir, sofa_get_data_fire,
%             intpolate_ir, ir_correct_distance

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
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


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 6;
narginchk(nargmin,nargmax)
if nargin==nargmax-1
    if isstruct(coordinate_system)
        conf = coordinate_system;
        coordinate_system = 'spherical';
    else
        conf = SFS_config;
    end
elseif nargin==nargmax-2
    coordinate_system = 'spherical';
    conf = SFS_config;
end
if length(head_orientation)==1
    head_orientation = [head_orientation 0];
end


%% ===== Computation ====================================================
%
% === SOFA loading ===
header = sofa_get_header(sofa);

% === Coordinate system conversion ===
% Convert everything to spherical coordinates
if strcmp('cartesian',coordinate_system)
    [xs(1),xs(2),xs(3)] = cart2sph(xs(1),xs(2),xs(3));
    [X(1),X(2),X(3)] = cart2sph(X(1),X(2),X(3));
elseif ~strcmp('spherical',coordinate_system)
    error('%s: unknown coordinate system type.',upper(mfilename));
end
% Store desired apparent position of source (in spherical coordinates)
xs = [correct_azimuth(xs(1)-X(1)-head_orientation(1)) ...
      correct_elevation(xs(2)-X(2)-head_orientation(2)) ...
      abs(X(3)-xs(3))];

% === Get Impulse Response ===
if strcmp('SimpleFreeFieldHRIR',header.GLOBAL_SOFAConventions)

    % http://www.sofaconventions.org/mediawiki/index.php/SimpleFreeFieldHRIR
    %
    % Returns a single impulse response for the specified position. The impulse
    % response is shifted in time and its amplitude is weighted according to the
    % desired distance. The desired direction is done by returning the nearest
    % neighbour or applying a linear interpolation.
    x0 = sofa_get_secondary_sources(header,'spherical');
    % Find nearest neighbours and interpolate if desired and needed
    [neighbours,idx] = findnearestneighbour(x0(:,1:2)',xs(1:2),3);
    ir = sofa_get_data_fir(sofa,idx);
    ir = ir_correct_distance(ir,x0(idx,3),xs(3),conf);
    ir = interpolate_ir(ir,neighbours,xs(1:2)',conf);

elseif strcmp('SingleRoomDRIR',header.GLOBAL_SOFAConventions)
    % FIXME: we should remove this altogether, as SingleRoomDRIR is original
    % considered to have only one source and use a microphone array:
    % http://www.sofaconventions.org/mediawiki/index.php/SingleRoomDRIR
    to_be_implemented;

elseif strcmp('MultiSpeakerBRIR',header.GLOBAL_SOFAConventions)
    %
    % http://www.sofaconventions.org/mediawiki/index.php/MultiSpeakerBRIR
    %
    % This looks for the nearest loudspeaker corresponding to the specified
    % source and listener position. For that loudspeaker the impulse reponse for
    % the specified head orientation is returned. If the head orientation could
    % not perfectly matched, an interpolation is applied as in the
    % SimpleFreeFieldHRIR case. If the specified head orientation is out of
    % bounds, the nearest head orientation is returned.
    %
    % Find nearest secondary source
    x0 = sofa_get_secondary_sources(header,'spherical');
    [~,idx_emitter] = findnearestneighbour(x0(:,1:3)',xs,1);
    % Find nearest head orientation in the horizontal plane
    [phi,theta] = sofa_get_head_orientations(header);
    [neighbours_head,idx_head] = ...
        findnearestneighbour([phi theta]',head_orientation,3);
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

else
    error('%s: %s convention is currently not supported.', ...
        upper(mfilename),header.GLOBAL_SOFAConventions);
end
