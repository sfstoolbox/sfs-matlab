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
%    x-axis                   
%       <--------------------------------------------------------
%                                   |
%          LS (Loudspeaker)         | 
%          |           [X0 Y0]      |
%          ^--^--^--^--^--^--^--^--^|-^--^                    
%                   (Array Center)  |
%                                   |
%                                   |
%                                   v y-axis
%
% see also: wfs_brs
%

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));

if ~isnumeric(L) || ~isscalar(L) || L<0
    error('%s: L has to be a positive scalar!');
end
if nargin<nargmax
    useconfig = true;
elseif ~isstruct(conf)
    error('%s: conf has to be a struct.',upper(mfilename));
else
    useconfig = false;
end


%% ===== Configuration ==================================================

% Load default configuration values
if(useconfig)
    conf = SFS_config;
end
% array type
array = conf.array;
% Center of the array
X0 = conf.X0;
Y0 = conf.Y0;
% Loudspeaker distance
LSdist = conf.LSdist;

%% ===== Calculation ====================================================

% Get the number of used loudspeaker
[nLS, L] = number_of_loudspeaker(L,conf);

if strcmp('linear',array)
    % === Linear loudspeaker array ===
    % Positions of the loudspeakers
    x0 = X0 + linspace(-L/2,L/2,nLS);
    y0 = Y0 * ones(1,nLS);
    % === Add jitter to the loudspeaker positions ===
    %x0(1) = X0-size/2;
    %x0(nLS) = X0+size/2;
    %for ii=2:nLS-1
    %    jitter = size/(4*nLS) * randn;
    %    x0(ii) = X0-size/2+(ii-1)*size/nLS + jitter;
    %end
    % Direction (orientation) of the loudspeaker
    phi = zeros(1,nLS);
elseif strcmp('circle',array)
    % === Circular loudspeaker array ===
    % Positions of the loudspeaker
    phi = linspace(-pi,(1-2/nLS)*pi,nLS);
    x0 = X0 + L/2*sin(-phi);
    y0 = Y0 + L/2*cos(-phi);
    % Direction of the loudspeakers
    phi = correct_angle(phi+pi);
elseif strcmp('box',array)
    % === Boxed loudspeaker array ===
    % Position and direction of the loudspeakers
    x0(1:nLS/4) = X0 + linspace(-L/2,L/2,nLS/4);
    y0(1:nLS/4) = Y0 + ones(1,nLS/4) * L/2 + LSdist;
    phi(1:nLS/4) = 0*ones(1,nLS/4);
    x0(2*nLS/4+1:3*nLS/4) = X0 + linspace(-L/2,L/2,nLS/4);
    y0(2*nLS/4+1:3*nLS/4) = Y0 - ones(1,nLS/4) * L/2 - LSdist;
    phi(2*nLS/4+1:3*nLS/4) = pi*ones(1,nLS/4);
    x0(nLS/4+1:2*nLS/4) = X0 + ones(1,nLS/4) * L/2 + LSdist;
    y0(nLS/4+1:2*nLS/4) = Y0 + linspace(-L/2,L/2,nLS/4);
    phi(nLS/4+1:2*nLS/4) = -pi/2*ones(1,nLS/4);
    x0(3*nLS/4+1:nLS) = X0 - ones(1,nLS/4) * L/2 - LSdist;
    y0(3*nLS/4+1:nLS) = Y0 + linspace(-L/2,L/2,nLS/4);
    phi(3*nLS/4+1:nLS) = pi/2*ones(1,nLS/4);
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
