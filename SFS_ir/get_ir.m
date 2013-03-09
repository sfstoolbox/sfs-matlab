function ir = get_ir(irs,xs,coordinate_system,conf)
%GET_IR returns a IR for the given apparent angle
%
%   Usage: ir = get_ir(irs,phi,[theta])
%
%   Input parameters:
%       irs                 - IR data set
%       xs                  - position of source in spherical (default) or
%                             cartesian coordinates
%       coordinate_system   - 'spherical' (default)
%                             'cartesian'
%       phi     - azimuth angle for the desired IR (rad)
%       theta   - elevation angle for the desired IR (rad)
%
%   Output parameters:
%       ir      - IR for the given source position (length of IR x 2)
%
%   GET_IR(irs,xs,coordinate_system) returns a single IR pair for the given
%   source position. The source position can be specified in spherical
%   coordinates in the form [phi theta r] in rad and meter. Alternativ you can
%   specify 'cartesian' as coordinate system and the position vector is then
%   [x y z9 in meter. If the desired angles are not present in the IR data set an
%   interpolation is applied to create the desired angles. The radius, that is
%   the distance of the source, is extrapolated by delaying and weighting with
%   1/r.
%
%   see also: read_irs, slice_irs, ir_intpol

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 4;
narginchk(nargmin,nargmax)
if nargin<nargmax
    conf = SFS_config;
end
if nargin<nargmax-1
    coordinate_system = 'spherical';
end
% Due to perfomance optimization the input parameter are not checked further


%% ===== Configuration ==================================================
% Precission of the wanted angle. If a IR within the given precission could be
% found no interpolation is applied.
prec = 0.1; % degree
c = conf.c;
fs = conf.fs;


%% ===== Computation ====================================================

% === Coordinate system conversion ===
if strcmp(coordinate_system,'cartesian') || strcmp(coordinate_system,'cart')
    [phi,theta,r] = cart2sph(xs(1),xs(2),xs(3));
else
    phi = xs(1);
    theta = xs(2);
    r = xs(3);
end
% === Check the given angles ===
% Ensure -pi <= phi < pi and -pi/2 <= theta <= pi/2
phi = correct_azimuth(phi);
theta = correct_elevation(theta);


% === IR interpolation ===
% Check if the IR dataset contains a measurement for the given angles
% phi and theta. If this is not the case, interpolate the dataset for the given
% angles.

% If azimuth and elevation could be found
idx = findrows(...
    roundto([irs.apparent_azimuth' irs.apparent_elevation'],prec),...
    roundto([phi,theta],prec));
if idx
    if length(idx)>1
        error(['%s: the irs data set has more than one entry corresponding ',...
               'an azimuth of %.3f deg and an elevation of %.3f deg.'],...
            upper(mfilename),degree(phi),degree(theta));
    end
    ir(:,1) = irs.left(:,idx);
    ir(:,2) = irs.right(:,idx);

% If only elevation angle could be found
elseif findrows(roundto(irs.apparent_elevation',prec), ...
                roundto(theta,prec))
    idx = findrows(roundto(irs.apparent_elevation',prec), ...
                   roundto(theta,prec));

    % === Interpolation of the azimuth ===
    % Get the IR set for the elevation theta
    irs = slice_irs(irs,idx);

    % Find the nearest value smaller than phi
    % Note: this requieres monotonic increasing values of phi in
    % azimuth(idx_theta)
    idx1 = find(irs.apparent_azimuth<phi,1,'last');
    if(isempty(idx1))
        % If no value is smaller than phi, use the largest value in
        % azimuth(idx_theta), because of the 0..2pi cicle
        idx1 = length(irs.apparent_azimuth(idx));
    end

    % Find the nearest value larger than phi
    idx2 = find(irs.apparent_azimuth>phi,1,'first');
    if(isempty(idx2))
        % If no value is greater than phi, use the smallest value in
        % azimuth(idx_theta), because of the 0..2pi cicle
        idx2 = 1;
    end

    if idx1==idx2
        error('%s: we have only one apparent_azimuth angle: %f.',...
            upper(mfilename),irs_theta.apparent_azimuth(idx1));
    end

    % Get the single IR corresponding to idx1
    ir1(:,1) = irs.left(:,idx1);
    ir1(:,2) = irs.right(:,idx1);
    % Get the single IR corresponding to idx2
    ir2(:,1) = irs.left(:,idx2);
    ir2(:,2) = irs.right(:,idx2);
    warning('SFS:irs_intpol',...
        ['doing IR interpolation with the angles beta1 = ',...
        '%.1f deg and beta2 = %.1f deg.'],...
        degree(irs.apparent_azimuth(idx1)),...
        degree(irs.apparent_azimuth(idx2)));
    % IR interpolation
    ir = intpol_ir(ir1,irs.apparent_azimuth(idx1),...
        ir2,irs.apparent_azimuth(idx2),phi);

% if only azimuth angle could be found
elseif findrows(roundto(irs.apparent_azimuth',prec),...
                roundto(phi,prec))
    idx = findrows(roundto(irs.apparent_azimuth',prec),...
                   roundto(phi,prec));

    % === Interpolation of the elevation ===
    % Get the IR set for the azimuth phi
    irs = slice_irs(irs,idx);

    % Find the nearest value smaller than theta
    % Note: this requieres monotonic increasing values of theta in
    % irs.apparent_theta(idx)
    idx1 = find(irs.apparent_elevation<theta,1,'last');
    if(isempty(idx1))
        error(['%s: there is no elevation avaialble which is smaller ',...
               'than your given theta value.'],upper(mfilename));
    end

    % Find the nearest value larger than theta
    idx2 = find(irs.apparent_elevation>theta,1,'first');
    if(isempty(idx2))
        error(['%s: there is no elevation avaialble which is greater ',...
               'than your given theta value.'],upper(mfilename));
    end

    if idx1==idx2
        error('%s: we have only one apparent_elevation angle: %f.',...
            upper(mfilename),irs_theta.apparent_azimuth(idx1));
    end

    % Get the single IR corresponding to idx1
    ir1(:,1) = irs.left(:,idx1);
    ir1(:,2) = irs.right(:,idx1);
    % Get the single IR corresponding to idx2
    ir2(:,1) = irs.left(:,idx2);
    ir2(:,2) = irs.right(:,idx2);
    % IR interpolation
    warning('SFS:irs_intpol',...
        ['doing IR interpolation with the angles beta1 = ',...
        '%.1f deg and beta2 = %.1f deg.'],...
        degree(irs.apparent_elevation(idx1)),...
        degree(irs.apparent_elevation(idx2)));
    ir = intpol_ir(ir1,irs.apparent_elevation(idx1),...
        ir2,irs.apparent_elevation(idx2),theta);

else
    error(['%s: at the moment interpolation for azimuth and elevation ',...
           'angles at the same time is currently not supported. ',...
           'Please choose an azimuth angle or an elevation angle, ',...
           'which is in the IR data set.'],upper(mfilename));
end


% === Distance handling ===
ir_distance = get_ir_distance(irs,xs);
% Define an offset to ensure r-ir_distance+offset > 0
% FIXME: is thsi really neccessary or should this be handled by the delayline()
% function?
offset = 3; % in m
% Time delay of the source (at the listener position)
delay = (r-ir_distance+offset)/c*fs; % in samples
% Amplitude weighting (point source model)
weight = 1/(4*pi*(r-ir_distance+offset));
% Apply delay and weighting
ir(:,1) = delayline(ir(:,1)',delay,weight,conf)';
ir(:,2) = delayline(ir(:,2)',delay,weight,conf)';


end % of main function

%% ===== Subfunctions ====================================================
% round the input matrix m to the given precission prec in degree
function m = roundto(m,prec)
    m = round(degree(m)/prec)*prec;
end
