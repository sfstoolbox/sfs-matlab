function ir = get_ir(sofa,X,head_orientation,xs,coordinate_system,conf)
%GET_IR returns an impulse response for the given apparent angle
%
%   Usage: ir = get_ir(sofa,X,head_orientation,xs,[coordinate_system],[conf])
%
%   Input parameters:
%       sofa    - impulse response data set (sofa struct/file)
%       X       - position of the listener, specified in the defined
%                 coordinate_system, see below
%       head_orientation  - orientation of the listener with phi / rad,
%                           delta / rad
%       xs      - position of the desired source, specified in the defined
%                 coordinate_system, see below
%       coordinate_system - coordinate system xs is specified in, avialable
%                           systems are:
%                             'spherical' - spherical system with phi / rad,
%                                           delta / rad, r / m  (default)
%                             'cartesian' - cartesian system with x,y,z / m
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       ir      - impulse response for the given position (length of IR x 2)
%
%   GET_IR(sofa,X,head_orientation,xs) returns a single impulse response from
%   the given SOFA file or struct. The impulse response is determined by the
%   position X and head orientation head_orientation of the listener, and the
%   position xs of the desired point source.
%   The desired distance between point source and listener is achieved by
%   delaying and weighting the impulse response. Distances larger than 10m are
%   ignored and set constantly to 10m. If the desired angles are not
%   present in the SOFA data set and conf.ir.useinterpolation is set to true
%   an interpolation is applied to create the impulse response.
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
%   FIXME: ask Piotr if we should now reference the AES standard together with
%   the web site?
%
%   see also: SOFAload, sofa_get_header, sofa_get_data, intpol_ir 

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


%% ===== Configuration ==================================================
% In subfunctions the following configurations are used, see below
% conf.ir.useinterpolation
% conf.ir.useoriglength
 

%% ===== Computation ====================================================
%
% === SOFA loading ===
header = sofa_get_header(sofa);

% === SOFA checking ===
% Conventions handled so far. For a list of SOFA conventions, see:
% http://www.sofaconventions.org/mediawiki/index.php/SOFA_conventions
SOFA_conventions = { ...
    'SimpleFreeFieldHRIR', ...
    'SingleRoomDRIR', ...
    };
if ~any(strcmp(header.GLOBAL_SOFAConventions,SOFA_conventions))
    error('%s: this SOFA Convention is currently not supported.');
end

% === Coordinate system conversion ===
if strcmp('cartesian',coordinate_system)
    [xs(1),xs(2),xs(3)] = cart2sph(xs(1),xs(2),xs(3));
    [X(1),X(2),X(3)] = cart2sph(X(1),X(2),X(3));
elseif ~strcmp('spherical',coordinate_system)
    error('%s: unknown coordinate system type.',upper(mfilename));
end
% Store desired apparent position of source
xs = [correct_azimuth(xs(1)-X(1)-head_orientation(1)) ...
      correct_elevation(xs(2)-X(2)-head_orientation(2)) ...
      abs(X(3)-xs(3))];

% === Get Impulse Response ===
if strcmp('SimpleFreeFieldHRIR',header.GLOBAL_SOFAConventions)
    % If we have free field HRTFs there is no difference between changing head
    % orientation or source poition. That means all we need is given by the
    % apparent source position stored in x0.
    ir = get_single_ir(sofa,xs,conf);

elseif strcmp('SingleRoomDRIR',header.GLOBAL_SOFAConventions)
    % For a given BRIR recording we are searching first for the correct head
    % orientation and look then for the source position
    %
    % Find indices for the given head orientation
    sofa_head_orientation = SOFAconvertCoordinates( ...
        header.ListenerView,header.ListenerView_Type,'spherical');
    tmp = vector_norm( ...
        bsxfun(@minus,rad(sofa_head_orientation(:,1:2)),head_orientation),2 ...
        );
    if min(tmp)>rad(1)  % Assuming resolution of 1 deg for head orientatons
        error('%s: no BRIR for head orientation (%i, %i) available', ...
            upper(mfilename), ...
            deg(head_orientation(1)), ...
            deg(head_orientation(2)));
    else
        idx = find(abs(tmp-min(tmp)<=eps));
        %sofa_head_orientation(idx,1:2)
        %[deg(x0(idx,1:2)) x0(idx,3)]
        %[deg(xs(1:2)); xs(3)]'
    end
    % Now, search for the correct source position
    ir = get_single_ir(sofa,x0,xs,conf);

end
end

function ir = correct_radius(ir,ir_distance,r,conf)
    % The basic idea is to handle the radius different then the directional
    % position. Between the directional positions an interpolation is applied if
    % desired. Different distances r are realized by weighting and adjusting the
    % time delay of the impulse response. As a starting value the actual
    % distance specified in the impulse response is used an adjusted regarding
    % the desired radius.
    % WARNING: in order to ensure correct extrapolation zeros corresponding to
    % the radius of the measured impulse response are padded at its beginning.
    % If you don't want this behavior you have to specify
    % conf.ir.useoriglength = true;
    % and ensure yourself that we have enough zeros at the beginning of each
    % impulse response and that possible different radii are represented with
    % the existing zeros.
    %
    % FIXME: check the performance of this function. The zeropadding at the
    % beginning of the ipulse response was done before only once for the whole
    % ir data set. Now, this is only possible if the complete sofa struct is
    % provided and not only the file. One could enhance the performance by doing
    % the zeropadding again for the whole database if the sofa struct is
    % provided. On the downside this would lead to a more complicated code base.
    %
    % Stop extrapolation for distances larger than 10m
    if any(ir_distance>10)
        ir_distance = min(ir_distance,10);
        warning(['%s: Your desired radius is larger than 10m, but we will ', ...
            'only extrapolate up to 10m. All larger radii will be set to ', ...
            '10m.'],upper(mfilename));
    end
    % Append zeros at the beginning of the impulse responses corresponding to
    % its maximum radius
    if ~conf.ir.useoriglength
        zero_padding = ceil(ir_distance/conf.c * conf.fs);
        if conf.N-zero_padding<128
            error(['%s: choose a larger conf.N value, because otherwise you ', ...
                'will have only %i samples of your original impulse response.'], ...
                upper(mfilename),conf.N-zero_padding);
        end
    else
        zero_padding = 0;
    end
    % Time delay of the source (at the listener position)
    delay = (r-ir_distance)/conf.c*conf.fs; % / samples
    % Amplitude weighting (point source model)
    % This gives weight=1 for r==ir_distance
    weight = ir_distance./r;
    if abs(delay)>size(ir,3)
        error(['%s: your impulse response is to short for a desired ', ...
            'delay of %.1f samples.'],upper(mfilename),delay);
    end
    % Apply delay and weighting
    ir = delayline(ir,[delay+zero_padding; delay+zero_padding], ...
        [weight; weight],conf);
end

function ir = get_ir_from_sofa(sofa,idx,x0,xs,conf)
    % Get the desired impulse response from the SOFA file/struct and returns it
    % with applied radius correction
    ir = sofa_get_data(sofa,idx);
    ir = correct_radius(ir,x0(idx,3),xs(3),conf);
end

function ir = get_single_ir(sofa,xs,conf)
    % GET_SINGLE_IR returns the desired impulse response from the given ones
    %
    % Input parameters:
    %   sofa    - SOFA file or struct
    %   xs      - position of the desired impulse response. xs is relative to
    %             [0 0 0] and given in spherical coordinates / [rad rad m].
    %   conf    - SFS configuration struct, specifying if interpolation between
    %             the given impulse responses shpuld be applied if the desired
    %             position is not given in the SOFA data set.
    %
    % Output parameters
    %   ir      - single impulse response

    % === Configuration ===
    useinterpolation = conf.ir.useinterpolation;
    % Precission of the wanted angle. If an impulse response within the given
    % precission could be found no interpolation is applied.
    prec = 0.001; % ~ 0.05 deg

    % === Computation ===
    % Calculate the apparent azimuth, elevation and distance
    x0 = SOFAcalculateAPV(sofa_get_header(sofa)); % / [deg deg m]
    % Convert angles to radian
    x0(:,1:2) = rad(x0(:,1:2));
    % Find the three nearest positions to the desired one (incorporating only angle
    % values and disregarding the radius)
    % FIXME: try to include radius here
    [neighbours,idx] = findnearestneighbour(x0(:,1:2)',xs(1:2),3);
    % Check if we have found directly the desired point or have to interpolate
    % bewteen different impulse responses
    if norm(neighbours(:,1)-xs(1:2)')<prec || ~useinterpolation
        disp('No interpolaton.');
        % Return the first nearest neighbour
        ir = get_ir_from_sofa(sofa,idx(1),x0,xs,conf);
    else
        % === IR interpolation ===
        % Check if we have to interpolate in one or two dimensions
        if norm(neighbours(1,1)-neighbours(1,2))<prec || ...
            norm(neighbours(2,1)-neighbours(2,2))<prec
            % --- 1D interpolation ---
            warning('SFS:irs_intpol',...
                ['doing 1D IR interpolation between (%.1f,%.1f) deg ', ...
                 'and (%.1f,%.1f) deg.'], ...
                deg(x0(idx(1),1)), deg(x0(idx(1),2)), ...
                deg(x0(idx(2),1)), deg(x0(idx(2),2)));
            ir = get_ir_from_sofa(sofa,idx(1:2),x0,xs,conf);
            ir = intpol_ir(ir,neighbours(:,1:2),xs(1:2)');
        else
            % --- 2D interpolation ---
            warning('SFS:irs_intpol3D',...
                ['doing 2D IR interpolation between (%.1f,%.1f) deg, ', ...
                 '(%.1f,%.1f) deg and (%.1f,%.1f) deg.'], ...
                deg(x0(idx(1),1)), deg(x0(idx(1),2)), ...
                deg(x0(idx(2),1)), deg(x0(idx(2),2)), ...
                deg(x0(idx(3),1)), deg(x0(idx(3),2)));
            ir = get_ir_from_sofa(sofa,idx(1:3),x0,r,conf);
            ir = intpot_ir(ir,neighbours(:,1:3),xs');
        end
    end
end

