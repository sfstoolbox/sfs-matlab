function [x0,y0,phi] = secondary_source_positions(L,conf)
%SECONDARY_SOURCE_POSITIONS Generates the loudspeaker positions and directions
%for the given loudspeaker array
%   Usage: [x0,y0,phi] = secondary_source_positions(L,conf)
%          [x0,y0,phi] = secondary_source_positions(L)
%
%   Input options:
%       L           - the size of the array (length for a linear array, diameter
%                     for a circle or a box
%       conf        - configuration struct
%
%   Output options:
%       x0, y0      - Position of the loudspeaker (m)
%       phi         - Directions of the loudspeakers (rad)
%
%   SECONDARY_SOURCES_POSITIONS(L) generates the loudspeaker
%   positions [x0,y0] and their directions phi for the given (conf.array)
%   loudspeaker array
%
%   Geometry (for the linear array):
%
%                                y-axis
%                                   ^
%                                   |
%                                   |
%                                   |
%          v--v--v--v--v--v--v--v--v|-v--v
%          |           [X0 Y0]      |
%    (Loudspeaker)  (Array center)  |
%                                   |
%       --------------------------------------------------------> x-axis
%
% see also: secondary_source_selection, tapwin
%

% AUTHOR: Sascha Spors, Hagen Wierstorf

% NOTE: If you wanted to add a new type of loudspeaker array, do it in a way,
% that the loudspeakers are ordered in a way, that one can go around for closed
% arrays. Otherwise the tapering window function will not work properly.


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargpositivescalar(L);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================

% array type
array = conf.array;
% Center of the array
X0 = conf.X0;
Y0 = conf.Y0;
% Loudspeaker distance
dx0 = conf.dx0;
% Given loudspeaker positions
x0 = conf.x0;
y0 = conf.y0;
phi = conf.phi;


%% ===== Calculation =====================================================

% Check if we have already predefined loudspeaker positions
if length(x0)>0
    isargvector(conf.x0,conf.y0,conf.phi);
    isargequallength(conf.x0,conf.y0,conf.phi);
    % If we have predefined loudspeaker positions leave this function
    return
end

% Get the number of used loudspeaker
[nls, L] = number_of_loudspeaker(L,conf);

if strcmp('linear',array)
    % === Linear loudspeaker array ===
    % Positions of the loudspeakers
    x0 = X0 + linspace(-L/2,L/2,nls);
    y0 = Y0 * ones(1,nls);
    % === Add jitter to the loudspeaker positions ===
    %x0(1) = X0-size/2;
    %x0(nls) = X0+size/2;
    %for ii=2:nls-1
    %    jitter = size/(4*nls) * randn;
    %    x0(ii) = X0-size/2+(ii-1)*size/nls + jitter;
    %end
    % Direction (orientation) of the loudspeaker
    % Default: pointing in y-direction
    phi = ones(1,nls)*pi/2;
elseif strcmp('circle',array)
    % === Circular loudspeaker array ===
    % Positions of the loudspeaker
    phi = linspace(0,(2-2/nls)*pi,nls); % 0..2pi
    [x0,y0] = sph2cart(phi,0,L/2);
    % Direction of the loudspeakers (it should pointing to the center, therefore
    % +pi)
    phi = correct_azimuth(phi+pi);
elseif strcmp('box',array)
    % === Boxed loudspeaker array ===
    % Position and direction of the loudspeakers
    % FIXME: nls/4 can be another value than an integer
    x0(1:nls/4) = X0 + linspace(-L/2,L/2,nls/4);
    y0(1:nls/4) = Y0 + ones(1,nls/4) * L/2 + dx0;
    phi(1:nls/4) = -pi/2 * ones(1,nls/4);
    x0(nls/4+1:2*nls/4) = X0 + ones(1,nls/4) * L/2 + dx0;
    y0(nls/4+1:2*nls/4) = Y0 + linspace(L/2,-L/2,nls/4);
    phi(nls/4+1:2*nls/4) = pi * ones(1,nls/4);
    x0(2*nls/4+1:3*nls/4) = X0 + linspace(L/2,-L/2,nls/4);
    y0(2*nls/4+1:3*nls/4) = Y0 - ones(1,nls/4) * L/2 - dx0;
    phi(2*nls/4+1:3*nls/4) = pi/2 * ones(1,nls/4);
    x0(3*nls/4+1:nls) = X0 - ones(1,nls/4) * L/2 - dx0;
    y0(3*nls/4+1:nls) = Y0 + linspace(-L/2,L/2,nls/4);
    phi(3*nls/4+1:nls) = 0 * ones(1,nls/4);
elseif strcmp('U',array)
    to_be_implemented(mfilename);
elseif strcmp('custom',array)
    x0 = conf.x0;
    y0 = conf.y0;
    phi = conf.phi;
    return;
else
    error('%s: %s is not a valid array type.',upper(mfilename),array);
end
