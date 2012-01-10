function irs = irs_with_particular_elevation(irs,delta)
%IRS_WITH_PARTICULAR_ELEVATION(irs,delta) returns an IRS set which only contains data for
%one elevation
%
%   Usage: irs = irs_with_particular_elevation(irs,delta)
%          irs = irs_with_particular_elevation(irs)
%
%   Input parameters:
%       irs     - IR data set
%       delta   - elevation angle for the desired IR (rad)
%                 default: 0
%
%   Output parameters:
%       irs      - IRS for the given elevation
%
%   IRS_WITH_PARTICULAR_ELEVATION(irs,delta) returns a IRS set for the given angle delta,
%   or by default the horizontal plane. Input should be an IRS-set with diffrent values
%   for the elevation
%
%   see also: slice_irs, new_irs

% AUTHOR: Lars-Erik Riechert, Hagen Wierstorf


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

% Finding the entries belonging to delta and slice the irs
idx = (( round(irs.apparent_elevation*10)==round(10*delta) ));
irs = slice_irs(irs,idx);
