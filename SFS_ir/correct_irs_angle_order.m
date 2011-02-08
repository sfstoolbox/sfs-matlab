function irs = correct_irs_angle_order(irs)
%CORRECT_IRS_ANGLE_ORDER reorders the angle entries of a irs to be increasing
%
%   Usage: irs = correct_irs_angle_order(irs)
%
%   Input options
%       irs - irs struct
%
%   Output options
%       irs - irs struct with corrected angle ordering
%
%   CORRECT_IRS_ANGLE_ORDER(irs) corrects the order of the azimuth and elevation
%   entries to be increasing over the whole range. This is needed for the
%   interpolation functions.
%
%   See also: ir_intpol

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
check_irs(irs);


%% ===== Computation =====================================================
% Sort azimuth angle and reorder the whole irs
[phi,idx] = sort(irs.apparent_azimuth);
irs.left = irs.left(:,idx);
irs.right = irs.right(:,idx);
irs.apparent_azimuth = irs.apparent_azimuth(idx);
irs.apparent_elevation = irs.apparent_elevation(idx);
if size(irs.head_azimuth)~=[1 1]
    irs.head_azimuth = irs.head_azimuth(idx);
end
if size(irs.head_elevation)~=[1 1]
    irs.head_elevation = irs.head_elevation(idx);
end
if size(irs.torso_azimuth)~=[1 1]
    irs.torso_azimuth = irs.torso_azimuth(idx);
end
if size(irs.torso_elevation)~=[1 1]
    irs.torso_elevation = irs.torso_elevation(idx);
end

%FIXME: try if it will work to reorder the whole irs afterwards regarding the
% apparent_elevation angle
