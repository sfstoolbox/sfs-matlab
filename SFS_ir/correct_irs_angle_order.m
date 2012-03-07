function irs = correct_irs_angle_order(irs,conf)
%CORRECT_IRS_ANGLE_ORDER reorders the angle entries of a irs to be increasing
%
%   Usage: irs = correct_irs_angle_order(irs,conf)
%          irs = correct_irs_angle_order(irs)
%
%   Input options
%       irs     - irs struct
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output options
%       irs     - irs struct with corrected angle ordering
%
%   CORRECT_IRS_ANGLE_ORDER(irs) corrects the order of the azimuth and elevation
%   entries to be increasing over the whole range. This is needed for the
%   interpolation functions.
%
%   See also: ir_intpol

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
if nargin==nargmax-1
    conf = SFS_config;
end
if conf.debug
    check_irs(irs);
end


%% ===== Computation =====================================================
% Sort azimuth angle and reorder the whole irs
[phi,idx] = sort(irs.apparent_azimuth);
irs.left = irs.left(:,idx);
irs.right = irs.right(:,idx);
irs.apparent_azimuth = irs.apparent_azimuth(idx);
irs.apparent_elevation = irs.apparent_elevation(idx);
if ~isequal(size(irs.head_azimuth),[1 1])
    irs.head_azimuth = irs.head_azimuth(idx);
end
if ~isequal(size(irs.head_elevation),[1 1])
    irs.head_elevation = irs.head_elevation(idx);
end
if ~isequal(size(irs.torso_azimuth),[1 1])
    irs.torso_azimuth = irs.torso_azimuth(idx);
end
if ~isequal(size(irs.torso_elevation),[1 1])
    irs.torso_elevation = irs.torso_elevation(idx);
end
if ~isequal(size(irs.distance),[1 1])
    irs.distance = irs.distance(idx);
end
if ~(isequal(size(irs.source_position),[3 1]) || ...
        isequal(size(irs.source_position),[1 3]))
    irs.source_position = irs.source_position(:,idx);
end

%FIXME: try if it will work to reorder the whole irs afterwards regarding the
% apparent_elevation angle
