function d = get_ir_distance(irs,phi,delta)
%GET_IR_DISTANCE returns the distance for the given apparent angle
%
%   Usage: ir = get_ir_distance(irs,phi,[delta])
%
%   Input parameters:
%       irs     - IR data set
%       phi     - azimuth angle for the desired IR (rad)
%       delta   - elevation angle for the desired IR (rad)
%
%   Output parameters:
%       d       - distace for the given angles
%
%   GET_IR_DISTANCE(irs,phi,delta) returns the distance for the given angles
%   phi and delta. If the desired angles are not present in the IR data set an
%   interpolation is applied to create the desired angles.
%
%   see also: get_ir, read_irs, slice_irs, ir_intpol

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
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
nargmax = 3;
narginchk(nargmin,nargmax)
if nargin==nargmax-1
    delta = 0;
end


%% ===== Computation ====================================================

% === Check the given angles ===
% Ensure -pi <= phi < pi and -pi/2 <= delta <= pi/2
phi = correct_azimuth(phi);
delta = correct_elevation(delta);

% Check if we have only one distance for the whole irs set and return that
% value
if isscalar(irs.distance)
    d = irs.distance;
    return;
end

% === IR interpolation ===
% Check if the IR dataset contains a measurement for the given angles
% phi and delta. If this is not the case, interpolate the dataset for the given
% angles.

% If we have found the both angles
% Precision of the conformance of the given angle and the desired one
prec = 0.1; % which is ca. 0.1 degree

% If azimuth and elevation could be found
idx = findrows(...
    roundto([irs.apparent_azimuth' irs.apparent_elevation'],prec),...
    roundto([phi,delta],prec));
if idx
    if length(idx)>1
        error(['%s: the irs data set has more than one entry corresponding ',...
               'an azimuth of %.3f and an elevation of %.3f.'],...
            upper(mfilename),degree(phi),degree(delta));
    end
    d = irs.distance(idx);

% If only the elevation angle is found
elseif findrows(roundto(irs.apparent_elevation',prec), ...
                roundto(delta,prec))
        idx = findrows(roundto(irs.apparent_elevation',prec), ...
                       roundto(delta,prec));

    % === Interpolation of the azimuth ===
    % Get the IR set for the elevation delta
    irs = slice_irs(irs,idx);

    % Find the nearest value smaller than phi
    % Note: this requieres monotonic increasing values of phi in
    % azimuth(idx_delta)
    idx1 = find(irs.apparent_azimuth<phi,1,'last');
    if(isempty(idx1))
        % If no value is smaller than phi, use the largest value in
        % azimuth(idx_delta), because of the 0..2pi cicle
        idx1 = length(irs.apparent_azimuth(idx));
    end

    % Find the nearest value larger than phi
    idx2 = find(irs.apparent_azimuth>phi,1,'first');
    if(isempty(idx2))
        % If no value is greater than phi, use the smallest value in
        % azimuth(idx_delta), because of the 0..2pi cicle
        idx2 = 1;
    end

    if idx1==idx2
        error('%s: we have only one apparent_azimuth angle: %f.',...
            upper(mfilename),irs_delta.apparent_azimuth(idx1));
    end

    % Get the single distance corresponding to idx1
    d1 = irs.distance(idx1);
    % Get the single distance corresponding to idx2
    d2 = irs.distance(idx2);
    warning('SFS:irs_intpol',...
        ['doing IR interpolation with the angles beta1 = ',...
        '%.1f deg and beta2 = %.1f deg.'],...
        degree(irs.apparent_azimuth(idx1)),...
        degree(irs.apparent_azimuth(idx2)));
    % distance interpolation
    d = d1 + (d2-d1) / ...
        (irs.apparent_azimuth(idx2)-irs.apparent_azimuth(idx1)) * ...
        (phi-irs.apparent_azimuth(idx1));

% If only the azimuth angle is found
elseif findrows(roundto(irs.apparent_azimuth',prec), ...
                roundto(phi,prec))
    idx = findrows(roundto(irs.apparent_azimuth',prec), ...
                   roundto(phi,prec));

    % === Interpolation of the elevation ===
    % Get the IR set for the azimuth phi
    irs = slice_irs(irs,idx);

    % Find the nearest value smaller than delta
    % Note: this requieres monotonic increasing values of delta in
    % irs.apparent_delta(idx)
    idx1 = find(irs.apparent_elevation<delta,1,'last');
    if(isempty(idx1))
        error(['%s: there is no elevation avaialble which is smaller ',...
               'than your given delta value.'],upper(mfilename));
    end

    % Find the nearest value larger than delta
    idx2 = find(irs.apparent_elevation>delta,1,'first');
    if(isempty(idx2))
        error(['%s: there is no elevation avaialble which is greater ',...
               'than your given delta value.'],upper(mfilename));
    end

    if idx1==idx2
        error('%s: we have only one apparent_elevation angle: %f.',...
            upper(mfilename),irs_delta.apparent_azimuth(idx1));
    end

    % Get the single distance corresponding to idx1
    d1 = irs.distance(idx1);
    % Get the single distance corresponding to idx2
    d2 = irs.distance(idx2);
    % IR interpolation
    warning('SFS:irs_intpol',...
        ['doing IR interpolation with the angles beta1 = ',...
        '%.1f deg and beta2 = %.1f deg.'],...
        degree(irs.apparent_elevation(idx1)),...
        degree(irs.apparent_elevation(idx2)));
    d = d1 + (d2-d1) / ...
        (irs.apparent_elevation(idx2)-irs.apparent_elevation(idx1)) * ...
        (phi-irs.apparent_elevation(idx1));

else
    error(['%s: at the moment interpolation for azimuth and elevation ',...
           'angles at the same time is currently not supported. ',...
           'Please choose an azimuth angle or an elevation angle, ',...
           'which is in the IR data set.'],upper(mfilename));
end
end % of main function

%% ===== Subfunctions ====================================================
% round the input matrix m to the given precission prec in degree
function m = roundto(m,prec)
    m = round(degree(m)/prec)*prec;
end
