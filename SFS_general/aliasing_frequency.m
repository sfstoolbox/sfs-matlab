function fal = aliasing_frequency(dx0,conf)
%ALIASING_FREQUENCY returns the aliasing frequency
%
%   Usage: fal = aliasing_frequency(dx0,conf)
%          fal = aliasing_frequency(dx0)
%
%   Input options:
%       dx0     - distance between adjacent secondary sources (m)
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output options:
%       fal     - aliasing frequency (Hz)
%
%   ALIASING_FREQUENCY(dx0) returns the aliasing frequency for the given
%   interspacing of secondary sources. The value is calculated after
%   spors2009.
%
%   S. Spors and J. Ahrens - Spatial sampling artifacts of wave field synthesis
%   for the reproduction of virtual point sources. 126th AES, May 2009.
%
%   see also: wave_field_mono_wfs_25d
%

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargpositivescalar(dx0);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
c = conf.c;


%% ===== Computation =====================================================
% FIXME: better calculation possible?
fal = c/(2*dx0);
