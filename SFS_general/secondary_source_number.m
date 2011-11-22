function [nls,L] = secondary_source_number(L,conf)
%SECONDARY_SOURCE_NUMBER calculate the number of loudspeaker for a linear WFS array
%   Usage: [nls,L] = secondary_source_number(L,conf)
%          [nls,L] = secondary_source_number(L)
%
%   Input parameters:
%       L       - length of the loudspeaker array (m)
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       nls     - number of needed loudspeaker
%       L       - real length of the loudspeaker array (corresponding to
%                 conf.dx0)
%
%   SECONDARAY_SOURCE_NUMBER(L,conf) calculates the number of needed loudspeaker for
%   the given array length L, using the config loudspeaker distance conf.dx0.
%   Also the real length L of such a loudspeaker array will returned. This is
%   neccessary because the user given length L is probably not compatible to the
%   given value conf.dx0.
%
%   see also: secondary_source_positions, secondary_source_selectio,
%   secondary_source_selection
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin)),
isargpositivescalar(L);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================

% Array type
array = conf.array;
% Loudspeaker distance
dx0 = conf.dx0;
% Predefined loudspeaker positions
x0 = conf.x0;
y0 = conf.y0;
phi = conf.phi;


%% ===== Calculation ====================================================
%
% Check if we have given loudspeaker positions
if length(x0>0)
    isargvector(conf.x0,conf.y0,conf.phi);
    isargequallength(conf.x0,conf.y0,conf.phi);
    nls = length(x0);
    L = L;
elseif strcmp('linear',array)
    % Number of loudspeaker
    nls = fix(L/dx0)+1;
    % Corresponding size of loudspeaker array
    L = (nls-1)*dx0;
elseif strcmp('circle',array)
    % L is the diameter!
    % Perimeter of the circle
    P = pi*L;
    % Number of loudspeakers
    %nls = fix(P/dx0)+1;
    nls = round(P/dx0);
    % Corresponding size of loudspeaker array
    L = (nls*dx0)/pi;
elseif strcmp('box',array)
    % FIXME: check what will happened with the loudspeakers on the edges!
    % Number of loudspeakers
    nls = 4*(fix(L/dx0)+1);
    % Corresponding size of loudspeaker array
    L = (nls/4-1)*dx0;
else
    error('%s: %s is a unknown array type.',upper(mfilename),array);
end
