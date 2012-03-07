function ir = get_ir(irs,phi,delta,conf)
%GET_IR returns a IR for the given apparent angle
%   Usage: ir = get_ir(irs,phi,delta,conf)
%          ir = get_ir(irs,phi,delta)
%          ir = get_ir(irs,phi)
%
%   Input parameters:
%       irs     - IR data set
%       phi     - azimuth angle for the desired IR (rad)
%       delta   - elevation angle for the desired IR (rad)
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       ir      - IR for the given angles (length of IR x 2)
%
%   GET_IR(irs,phi,delta) returns a single IR set for the given angles phi and
%   delta. If the desired angles are not present in the IR data set an 
%   interpolation is applied to create the desired angles.
%
%   see also: read_irs, slice_irs, ir_intpol
%

% AUTHOR: Sascha Spors, Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 4;
error(nargchk(nargmin,nargmax,nargin))
if nargin==nargmax-1
    conf = SFS_config;
elseif nargin==nargmax-2
    delta = 0;
    conf = SFS_config;
end
if conf.debug
    check_irs(irs);
    isargscalar(phi,delta);
end


%% ===== Computation ====================================================

% === Check the given angles ===
% Ensure -pi <= phi < pi and -pi/2 <= delta <= pi/2
phi = correct_azimuth(phi,conf);
delta = correct_elevation(delta,conf);

% === IR interpolation ===
% Check if the IR dataset contains a measurement for the given angles
% phi and delta. If this is not the case, interpolate the dataset for the given
% angles.

% If we have found the both angles
% Precision of the conformance of the given angle and the desired one
prec = 10; % which is 0.1 degree
if findrows(...
    round(degree(prec*[irs.apparent_azimuth' irs.apparent_elevation'], ...
        conf)),...
    round(degree(prec*[phi,delta],conf)))

    idx = findrows(...
        round(degree(prec*[irs.apparent_azimuth' irs.apparent_elevation'], ...
            conf)),...
        round(degree(prec*[phi,delta],conf)));
    if length(idx)>1
        error(['%s: the irs data set has more than one entry corresponding ',...
               'an azimuth of %f and an elevation of %f.'],...
            upper(mfilename),phi,delta);
    end
    ir(:,1) = irs.left(:,idx);
    ir(:,2) = irs.right(:,idx);

elseif findrows(irs.apparent_elevation',delta)
    idx = findrows(irs.apparent_elevation',delta,conf);
    % === Interpolation of the azimuth ===
    % Get the IR set for the elevation delta
    irs = slice_irs(irs,idx,conf);

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

    % Get the single IR corresponding to idx1
    ir1(:,1) = irs.left(:,idx1);
    ir1(:,2) = irs.right(:,idx1);
    % Get the single IR corresponding to idx2
    ir2(:,1) = irs.left(:,idx2);
    ir2(:,2) = irs.right(:,idx2);
    warning('SFS:irs_intpol',...
        ['doing IR interpolation with the angles beta1 = ',...
        '%.1f deg and beta2 = %.1f deg.'],...
        degree(irs.apparent_azimuth(idx1),conf),...
        degree(irs.apparent_azimuth(idx2),conf));
    % IR interpolation
    ir = intpol_ir(ir1,irs.apparent_azimuth(idx1),...
        ir2,irs.apparent_azimuth(idx2),phi,conf);

elseif findrows(irs.apparent_azimuth',phi)
    idx = findrows(irs.apparent_azimuth',phi,conf);
    % === Interpolation of the elevation ===
    % Get the IR set for the azimuth phi
    irs = slice_irs(irs,idx,conf);

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
        degree(irs.apparent_elevation(idx1),conf),...
        degree(irs.apparent_elevation(idx2),conf));
    ir = intpol_ir(ir1,irs.apparent_elevation(idx1),...
        ir2,irs.apparent_elevation(idx2),delta,conf);

else
    error(['%s: at the moment interpolation for azimuth and elevation ',...
           'angles at the same time is currently not supported. ',...
           'Please choose an azimuth angle or an elevation angle, ',...
           'which is in the IR data set.'],upper(mfilename));
end
