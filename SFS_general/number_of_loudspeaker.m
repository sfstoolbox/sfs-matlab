function [nLS L] = number_of_loudspeaker(L,conf)
%NUMBER_OF_LOUDSPEAKER calculate the number of loudspeaker for a linear WFS array
%   Usage: nLS = number_of_loudspeaker(L,conf)
%          nLS = number_of_loudspeaker(L)
%
%   Input parameters:
%       L       - length of the loudspeaker array (m)
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       nLS     - number of needed loudspeaker
%       L       - real length of the loudspeaker array (correspnding to
%                 conf.LSdist)
%
%   NUMBER_OF_LOUDSPEAKER(L,conf) calculates the number of needed loudspeaker for
%   the given array length L, using the config loudspeaker distance conf.LSdist.
%   Also the real length L of such a loudspeaker array will returned. This is
%   neccessary because the user given length L is probably not compatible to the
%   given value conf.LSdist.
%
%   see also: secondary_source_positions
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin)),
isargpositivescalar({L},{'L'});

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct({conf},{'conf'});
end


%% ===== Configuration ==================================================

% Array type
array = conf.array;

% Loudspeaker distance
dLS = conf.LSdist;

%% ===== Calculation ====================================================
%
if strcmp('linear',array)
    % Number of loudspeaker
    nLS = fix(L/dLS)+1;
    % Corresponding size of loudspeaker array
    L = (nLS-1)*dLS;
elseif strcmp('circle',array)
    % L is the radius!
    % Perimeter of the circle
    P = pi*L;
    % Number of loudspeakers
    nLS = fix(P/dLS)+1;
    % Corresponding size of loudspeaker array
    L = ((nLS-1)*dLS)/pi;
elseif strcmp('box',array)
    % FIXME: check what will happened with the loudspeakers on the edges!
    % Number of loudspeakers
    nLS = 4*(fix(L/dLS)+1);
    % Corresponding size of loudspeaker array
    L = (nLS/4-1)*dLS;
else
    error('%s: %s is a unknown array type.',upper(mfilename),array);
end
