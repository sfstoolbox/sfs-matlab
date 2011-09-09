function d = get_ir_distance(irs,phi,delta)
%GET_IR_DISTANCE returns the distance for the given apparent angle
%   Usage: ir = get_ir_distance(irs,phi,delta)
%          ir = get_ir_distance(irs,phi)
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
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin))
check_irs(irs);
isargscalar(phi);
if nargin==nargmax
    isargscalar(delta);
else
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
prec = 1000; % which is ca. 0.1 degree
if findrows(round(prec*[irs.apparent_azimuth' irs.apparent_elevation']),...
        round(prec*[phi,delta]))
    idx = findrows(round(prec*[irs.apparent_azimuth' irs.apparent_elevation']),...
        round(prec*[phi,delta]));
    if length(idx)>1
        error(['%s: the irs data set has more than one entry corresponding ',...
               'an azimuth of %f and an elevation of %f.'],...
            upper(mfilename),phi,delta);
    end
    d = irs.distance(idx);

elseif findrows(irs.apparent_elevation',delta)
    idx = findrows(irs.apparent_elevation',delta);
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

elseif findrows(irs.apparent_azimuth',phi)
    idx = findrows(irs.apparent_azimuth',phi);
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
