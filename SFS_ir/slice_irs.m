function irspart = slice_irs(irs,idx)
%SLICE_IRS returns a part of an IRs set given by idx
%   Usage: irspart = slice_irs(irs,idx)
%
%   Input parameters:
%       irs     - IR data set
%       idx     - idx to slice out of the IR set
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


%% ===== Checking of input  parameters ==================================

if nargchk(2,2,nargin)
    error('Wrong number of args. Usage: irspart = slice_ir(irs,idx)');
end

if ~isstruct(irs)
    error('%s: irs has to be a struct!',upper(mfilename));
end
if ~isvector(idx)
    error('%s: idx has to be an index vector.',upper(mfilename));
end
check_irs(irs);

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
