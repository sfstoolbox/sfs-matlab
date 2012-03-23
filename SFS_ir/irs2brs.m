function brs = irs2brs(irs)
%IRS2BRS converts a irs data set to an brs set suitable for the SSR
%
%   Usage: brs = irs2brs(irs)
%
%   Input parameters:
%       irs     - irs data set
%
%   Output parameters:
%       brs     - brs data set
%
%   IRS2BRS(irs) converts a irs data set into a brs set suitable for the
%   SoundScape Renderer. The brs data set is a matrix containing the
%   channels for all directions.
%
%   see also: 

% AUTHOR: Hagen Wierstorf
% $LastChangedDate:$
% $LastChangedRevision:$
% $LastChangedBy:$


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin))
if nargin==nargmax-1
    conf = SFS_config;
end
check_irs(irs);
isargstruct(conf);


%% ===== Main ===========================================================

% Check if only one elevation angle is given
if length(unique(irs.apparent_elevation))~=1
    error(['%s: Your irs set has different elevation angles, which is',...
        ' not supported by the SoundScape Renderer.'],upper(mfilename));
end

% TODO: check the order of angles
%       I think the user have to check this by itself. Because the user
%       could also be interested in a particular angle order, they can
%       create by themself

for ii = 1:length(irs.apparent_azimuth)
    brs(:,(ii-1)*2+1:ii*2) = [irs.left(:,ii) irs.right(:,ii)];
end


        