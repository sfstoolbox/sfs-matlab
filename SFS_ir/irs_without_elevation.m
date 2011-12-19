function irs = irs_without_elevation(irs,delta)
%IRS_WITHOUT_ELEVATION(irs,delta) returns an IRS set which only contains data for one elevation
%   Usage: irs = irs_without_elevation(irs,delta)
%          irs = irs_without_elevation(irs)
%
%   Input parameters:
%       irs     - IR data set
%       delta   - elevation angle for the desired IR (rad)
%
%   Output parameters:
%       irs      - IRS for the given elevation 
%
%   IRS_WITHOUT_ELEVATION(irs,delta) returns a IRS set for the given angle delta, or by default
%	the horizontal plane. Input should be an IRS-set with diffrent values for the elevation
%
% AUTHOR: Lars-Erik Riechert

%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin))
check_irs(irs);
if nargin==nargmax
    isargscalar(delta);
else
    delta = 0;
end


%% ===== Computation ====================================================

%finding the entries belonging to delta
indexes = find(round(irs.apparent_elevation*10)) == round(10*delta);

irs2 = new_irs();

irs2.description = irs.description;
irs2.head = irs.head;
irs2.ears = irs.ears;
irs2.room = irs.room;
irs2.room_corners = irs.room_corners;
irs2.source = irs.source;
irs2.distance = irs.distance;
irs2.fs = irs.fs;

for i = 1:length(indexes)
	irs.apparent_elevation(:) = 
	irs.left(:,:) = irs.left(:,indexes)
irs.right(:,:) = irs.right(:,indexes)
end



check_irs(irs);