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
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
if ~isstruct(irs)
    error('%s: irs has to be a struct!',upper(mfilename));
end
% Check if the given irs is in the right format
if ~isvector(irs.azimuth) || ~isvector(irs.elevation)
    error(['%s: the given irs is not in the right format. ', ...
           'irs.azimuth and irs.elevation has to be avaiable.'], ...
        upper(mfilename));
end


%% ===== Computation =====================================================
% Reorder azimuth angle
phi = irs.azimuth;
[phi,idx] = sort(phi);
irs.azimuth = irs.azimuth(idx);
irs.elevation = irs.elevation(idx);
irs.left = irs.left(:,idx);
irs.right = irs.right(:,idx);

% FIXME: Reorder elevation angles for a given azimuth angle
