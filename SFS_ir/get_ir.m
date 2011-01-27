function ir = get_ir(irs,phi)
%GET_IR returns a HRIR/BRIR for the given angle
%   Usage: ir = get_ir(irs,phi)
%
%   Input parameters:
%       irs     - HRIR/BRIR data set for the second sources
%       phi     - angle for the desired IR (rad)
%
%   Output parameters:
%       ir      - HRIR/BRIR for the given angle phi (length of IR x 2)
%
%   GET_IR(irs,phi) returns a single IR set for the given angle phi. If the
%   desired angle is not present in the IR data set an interpolation is applied
%   to create the desired angle.
%
%   see also: read_irs, slice_irs, ir_intpol
%

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================

error(nargchk(2,2,nargin))

if ~isstruct(irs)
    error('%s: hrirs has to be a struct!',upper(mfilename));
end
% Check if we have only an azimuth value or an elevation as well
if length(phi)==2
    delta = phi(2);
    phi = phi(1);
else
    delta = 0;
end
if ~isnumeric(phi) || ~isscalar(phi)
    error('%s: phi has to be a scalar!',upper(mfilename));
end
if ~isnumeric(delta) || ~isscalar(delta)
    error('%s: delta has to be a scalar!',upper(mfilename));
end


%% ===== Computation ====================================================

% === Get the angles from the IR set ===
% Angles of the IR set
azimuth = irs.angle(1,:);
elevation = irs.angle(2,:);

% === Check the given angles ===
% FIXME: this works at the moment only for azimuth angles
% Ensure -pi <= phi < pi and -pi/2 <= delta < pi/2
%[phi,delta] = correct_angle([phi,delta]);
phi = correct_angle(phi);


% === IR interpolation ===
% Check if the IR dataset contains a measurement for the given angles
% phi and delta. If this is not the case, interpolate the dataset for the given
% angles.
idx_phi = (( azimuth==phi ));
idx_delta = (( elevation==delta ));

% Check if we have found the desired angles
if strmatch([phi,delta], irs.angle')

    % If no interpolation is needed use the desired IR
    idx = strmatch([phi,delta], irs.angle');
    ir(:,1) = irs.left(:,idx);
    ir(:,2) = irs.right(:,idx);

elseif ~sum(idx_phi) && sum(idx_delta)

    % === Interpolation of the azimuth ===

    % Get the IR set for the elevation delta
    irs_delta = slice_irs(irs,idx_delta);

    % Find the nearest value smaller than phi
    % Note: this requieres monotonic increasing values of phi in
    % azimuth(idx_delta)
    idx1 = find(azimuth(idx_delta) < phi,1,'last');
    if(isempty(idx1))
        % If no value is smaller than phi, use the largest value in
        % azimuth(idx_delta), because of the 0..2pi cicle
        idx1 = length(azimuth(idx_delta));
    end

    % Get the single IR corresponding to idx1
    ir1 = slice_irs(irs_delta,idx1);


    % Find the nearest value larger than phi
    idx2 = find(azimuth(idx_delta) > phi,1,'first');
    if(isempty(idx2))
        % If no value is greater than phi, use the smallest value in
        % azimuth(idx_delta), because of the 0..2pi cicle
        idx2 = 1;
    end

    % Get the single IR corresponding to idx2
    ir2 = slice_irs(irs_delta,idx2);

    % IR interpolation
    ir = ir_intpol(ir1,ir2,phi);

elseif sum(idx_phi) && ~sum(idx_delta)

    % TODO: this can not work probably, because of the missing 0..2pi
    % circle we have for the elevations!

    % === Interpolation of the elevation ===

    % Get the IR set for the elevation delta
    irs_phi = slice_irs(irs,idx_phi);

    % Find the nearest value smaller than delta
    % Note: this requieres monotonic increasing values of delta in
    % elevation(idx_phi)
    idx1 = find(elevation(idx_phi) < delta,1,'last');
    if(isempty(idx1))
        % If no value is smaller than delta, use the largest value in
        % elevation(idx_phi), because of the 0..2pi cicle
        idx1 = length(elevation(idx_phi));
    end

    % Get the single IR corresponding to idx1
    ir1 = slice_irs(irs_phi,idx1);


    % Find the nearest value larger than delta
    idx2 = find(elevation(idx_phi) > delta,1,'first');
    if(isempty(idx2))
        % If no value is greater than delta, use the smallest value in
        % elevation(idx_phi), because of the 0..2pi cicle
        idx2 = 1;
    end

    % Get the single IR corresponding to idx2
    ir2 = slice_irs(irs_phi,idx2);

    % IR interpolation
    ir = ir_intpol(ir1,ir2,delta,conf);

else

    error(['%s: The elevation angle or the azimuth angle has to be avaiable',...
           'in the IR data set in order to apply an interpolation!.'],...
           upper(mfilename));

end
