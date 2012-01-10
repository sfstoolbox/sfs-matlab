function irs = irs_with_particular_azimuth(irs,phi)
%IRS_WITH_PARTICULAR_ELEVATION(irs,phi) returns an IRS set which only contains data for
%one azimuth
%
%   Usage: irs = irs_with_particular_azimuth(irs,phi)
%          irs = irs_with_particular_azimuth(irs)
%
%   Input parameters:
%       irs     - IR data set
%       phi     - azimuth angle for the desired IR (rad)
%                 default: 0
%
%   Output parameters:
%       irs      - IRS for the given azimuth
%
%   IRS_WITH_PARTICULAR_ELEVATION(irs,phi) returns a IRS set for the given angle phi,
%   or by default the medial plane. Input should be an IRS-set with diffrent values
%   for the azimuth
%
%   see also: slice_irs, new_irs

% AUTHOR: Lars-Erik Riechert, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin))
check_irs(irs);
if nargin==nargmax
    isargscalar(phi);
else
    phi = 0;
end


%% ===== Computation ====================================================

% Finding the entries belonging to phi and slice the irs
idx = (( round(irs.apparent_azimuth*10)==round(10*phi) ));
irs = slice_irs(irs,idx);
