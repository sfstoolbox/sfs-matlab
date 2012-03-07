function irspart = slice_irs(irs,idx,conf)
%SLICE_IRS returns a part of an IRs set given by idx
%   Usage: irspart = slice_irs(irs,idx,conf)
%          irspart = slice_irs(irs,idx)
%
%   Input parameters:
%       irs     - IR data set
%       idx     - idx to slice out of the IR set
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       irspart - HRIR/BRIR containing only the part of the original IR set 
%                 given by idx
%
%   SLICE_IRS(irs,idx) returns a part of the IR set irs given by idx. The new 
%   part is a full IR set containing all the necessary struct elements.
%
%   see also: get_ir, read_irs
%

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));
if nargin==nargmax-1
    conf = SFS_config;
end
if conf.debug
    isargvector(idx);
    check_irs(irs);
end


%% ===== Slicing the IR set ==============================================
irspart = irs;
irspart.left = irs.left(:,idx);
irspart.right = irs.right(:,idx);
irspart.apparent_azimuth = irs.apparent_azimuth(idx);
irspart.apparent_elevation = irs.apparent_elevation(idx);
if size(irs.head_azimuth)~=[1 1]
    irspart.head_azimuth = irs.head_azimuth(idx);
end
if size(irs.head_elevation)~=[1 1]
    irspart.head_elevation = irs.head_elevation(idx);
end
if size(irs.torso_azimuth)~=[1 1]
    irspart.torso_azimuth = irs.torso_azimuth(idx);
end
if size(irs.torso_elevation)~=[1 1]
    irspart.torso_elevation = irs.torso_elevation(idx);
end
if size(irs.distance)~=[1 1]
    irspart.distance = irs.distance(idx);
end
if size(irs.source_position)~=[3 1]
    irspart.source_position = irs.source_position(idx);
end
if size(irs.source_reference)~=[3 1]
    irspart.source_reference = irs.source_reference(idx);
end
