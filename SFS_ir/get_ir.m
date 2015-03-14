function ir = get_ir(sofa,X,head_orientation,xs,coordinate_system,conf)
%GET_IR returns an impulse response for the given apparent angle
%
%   Usage: ir = get_ir(sofa,X,head_orientation,xs,[coordinate_system],[conf])
%
%   Input parameters:
%       sofa    - impulse response data set (sofa struct/file)
%       X       - position of the listener, specified in the defined
%                 coordinate_system, see below
%       xs      - position of the desired source, specified in the defined
%                 coordinate_system, see below
%       head_orientation  - orientation of the listener with phi / rad,
%                           delta / rad
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
%   GET_IR(sofa,phi,delta,r,X0) returns a single IR set for the given angles 
%   phi and delta. If the desired angles are not present in the IR data set 
%   an interpolation is applied to create the desired angles.
%   The desired radius is achieved by delaying and weighting the impulse
%   response.
%   For a description of the SOFA file format see: http://sofaconventions.org
%   FIXME: ask Piotr if we should now reference the AES standard together with
%   the web site?
%
%   see also: read_irs, slice_irs, intpol_ir 

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
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
useinterpolation = conf.ir.useinterpolation;
% Precission of the wanted angle. If an impulse response within the given
% precission could be found no interpolation is applied.
prec = 0.001; % ~ 0.05 deg


%% ===== Computation ====================================================

% === SOFA loading ===
% Get only the metadata of the SOFA data set
if ~isstruct(sofa) && exist(sofa,'file')
    % Get metadata
    header = SOFAload(sofa,'nodata');
elseif isfield(sofa.Data,'IR')
    header = sofa;
    header.Data = rmfield(sofa.Data,'IR');
else
    error('%s: sofa has to be a file or a SOFA struct.',upper(mfilename));
end

% === SOFA checking ===
% Conventions handled so far (see: http://...FIXME)
SOFA_conventions = { ...
    'SimpleFreeFieldHRIR', ...
    'SingleRoomDRIR', ...
    };
if ~any(strcmp(header.GLOBAL_SOFAConventions,SOFA_conventions))
    error('%s: wrong SOFA Convention.');
end

% === Internal variables ===
% Calculate the apparent azimuth, elevation and distance
x0 = SOFAcalculateAPV(header); % / [deg deg m]
% Convert angles to radian
x0(:,1:2) = rad(x0(:,1:2));


% === Coordinate system conversion ===
if strcmp('cartesian',coordinate_system)
    [xs(1),xs(2),xs(3)] = cart2sph(xs(1),xs(2),xs(3));
    [X(1),X(2),X(3)] = cart2sph(X(1),X(2),X(3));
elseif ~strcmp('spherical',coordinate_system)
    error('%s: unknown coordinate system type.',upper(mfilename));
end
% Store desired relative position of source
r = abs(X(3)-xs(3));
xs = [correct_azimuth(xs(1)-X(1)-head_orientation(1)) ...
      correct_elevation(xs(2)-X(2)-head_orientation(2))]';

if strcmp('SimpleFreeFieldHRIR',header.GLOBAL_SOFAConventions)
    disp('SimpleFreeFieldHRIR')
    % If we have free field HRTFs there is no difference between changing head
    % orientation or source poition. That means all we need is given by the
    % apparent source position stored in x0.
    

    % Find the three nearest positions to the desired one (incorporating only angle
    % values and disregarding the radius)
    % FIXME: try to include radius here
    [neighbours,idx] = findnearestneighbour(x0(:,1:2)',xs,3);

    % Check if we have found directly the desired point or have to interpolate
    % bewteen different impulse responses
    if norm(neighbours(:,1)-xs)<prec || ~useinterpolation
        % Return the first nearest neighbour
        ir = get_sofa_data(sofa,idx(1));
        ir = correct_radius(ir,x0(idx(1),3),r,conf);
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
            ir1 = get_sofa_data(sofa,idx(1));
            ir2 = get_sofa_data(sofa,idx(2));
            ir3 = get_sofa_data(sofa,idx(3));
            ir1 = correct_radius(ir1,x0(idx(1),3),r,conf);
            ir2 = correct_radius(ir2,x0(idx(2),3),r,conf);
            ir = intpol_ir(ir1,ir2,neighbours(:,1:2),xs);
        else
            % --- 2D interpolation ---
            warning('SFS:irs_intpol3D',...
                ['doing 2D IR interpolation between (%.1f,%.1f) deg, ', ...
                 '(%.1f,%.1f) deg and (%.1f,%.1f) deg.'], ...
                deg(x0(idx(1),1)), deg(x0(idx(1),2)), ...
                deg(x0(idx(2),1)), deg(x0(idx(2),2)), ...
                deg(x0(idx(3),1)), deg(x0(idx(3),2)));
            ir1 = get_sofa_data(sofa,idx(1));
            ir2 = get_sofa_data(sofa,idx(2));
            ir3 = get_sofa_data(sofa,idx(3));
            ir1 = correct_radius(ir1,x0(idx(1),3),r,conf);
            ir2 = correct_radius(ir2,x0(idx(2),3),r,conf);
            ir3 = correct_radius(ir3,x0(idx(3),3),r,conf);
            ir = intpot_ir(ir1,ir2,ir3,[neighbours; 1 1 1],[xs;1]);
        end
    end

elseif strcmp('SingleRoomDRIR',header.GLOBAL_SOFAConventions)
    disp('SingleRoomDRIR')
    % For a given BRIR recording we are searching first for the correct head
    % orientation and look then for the source position

    % Find indices for the given head orientation
    sofa_head_orientation = SOFAconvertCoordinates( ...
        header.ListenerView,header.ListenerView_Type,'spherical');
    tmp = vector_norm( ...
        bsxfun(@minus,rad(sofa_head_orientation(:,1:2)),head_orientation),2);
    if min(tmp)>rad(2)
        error('%s: no BRIR for head orientation (%i, %i) available', ...
            upper(mfilename), ...
            deg(head_orientation(1)), ...
            deg(head_orientation(2)));
    else
        idx = ( tmp==min(tmp) );
        x0 = x0(idx,:);
    end

    % Now, search for the correct source position
    [~,idx] = findnearestneighbour(x0(:,1:2)',xs,1);
    ir = get_sofa_data(sofa,idx(1));

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
    if ir_distance>10
        ir_distance = 10;
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
    % Delay only if we have an delay other than 0
    if abs(r-ir_distance)>0.0001 % ~0.01 samples
        % Time delay of the source (at the listener position)
        delay = (r-ir_distance)/conf.c*conf.fs; % / samples
        % Amplitude weighting (point source model)
        % This gives weight=1 for r==ir_distance
        weight = 1/(r-ir_distance+20) * 20;
        if abs(delay)>size(ir,1)
            error(['%s: your impulse response is to short for a desired ', ...
                'delay of %i samples.'],upper(mfilename),delay);
        end
    else
        delay = 0;
        weight = 1;
    end
    % Apply delay and weighting
    ir = delayline(ir,[delay+zero_padding; delay+zero_padding], ...
        [weight; weight],conf);
end

function ir = get_sofa_data(sofa,idx)
    % Check if we have a file or struct
    if isstruct(sofa)
        ir = sofa.Data.IR(idx,:,:);
        ir = squeeze(ir)';
    else
        file = sofa;
        % FIXME: this is currently not working under Octave, because it will
        % always return [1 idx] data
        sofa = SOFAload(file,[idx idx]);
        ir = squeeze(sofa.Data.IR)';
    end
end
