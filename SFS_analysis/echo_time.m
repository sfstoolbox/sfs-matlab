function t = echo_time(X,Y,xs,ys,L,conf)
%ECHO_TIME time of occurence of echo for a linear WFS array
%   Usage: echo_time(X,Y,xs,ys,L,conf)
%          echo_time(X,Y,xs,ys,L)
%
%   Input parameters:
%       X,Y     - listener position
%       xs,ys   - position of the virtual source
%       L       - length of the linear loudspeaker array
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   ECHO_TIME(X,Y,xs,ys,L) calculates the time of occuring of the first
%   echo (focused sources) resp. the last echo (virtual point source) at the
%   listener position X,Y for a given point source location xs,ys and a 
%   linear WFS loudspeaker array with a length of L.
%
%   see also: plot_echo_times
%

% AUTHOR: Hagen Wierstorf

%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 6;
error(nargchk(nargmin,nargmax,nargin));

if ~isnumeric(X) || ~isscalar(X)
    error('%s: X has to be a scalar!',upper(mfilename));
end

if ~isnumeric(Y) || ~isscalar(Y)
    error('%s: Y has to be a scalar!',upper(mfilename));
end

if ~isnumeric(xs) || ~isscalar(xs)
    error('%s: xs has to be a scalar!',upper(mfilename));
end

if ~isnumeric(ys) || ~isscalar(ys)
    error('%s: ys has to be a scalar!',upper(mfilename));
end

if ~isnumeric(L) || ~isscalar(L) || L<0
    error('%s: L has to be a positive scalar!',upper(mfilename));
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

% Loudspeaker distance
LSdist = conf.LSdist;
% Array center position
X0 = conf.X0;
Y0 = conf.Y0;

% Speed of sound
c = conf.c;

% Bandwidth of Dirac pulse
%f = fix(fs/2);
%k = 2*pi*f/c;


%% ===== Variables ======================================================
% Number of loudspeaker (round towards plus infinity)
nLS = ceil(L/LSdist);

% Loudspeaker positions
[LSpos,LSdir] = LSpos_linear(X0,Y0,(nLS-1)*LSdist,nLS);
x0 = LSpos(1,:);
y0 = LSpos(2,:);


%% ===== Calculate a time axis ==========================================

% Geometry
%          [x0 y0]                [X0 Y0]
% x-axis <-^--^--^--^--^--^--^--^--^-|-^--^--^--^--^--^--^--^--^--
%             |                      |
%         R2 |  |                    |
%           |     | R                |      
%          x        |                |
%       [xs,ys]       |              |
%                       O            |
%                      [X,Y]         |
%                                    v
%                                  y-axis
%
% Distance between secondary source and listener: R
% Distance between secondary source and virtual source: R2

% Calculate arrival times at the virtual source position of the waves 
% emitted by the secondary sources
t2 = zeros(nLS,1);
t = zeros(nLS,1);
for i = 1:nLS
    % Distance between secondary sources and virtual source
    R2 = norm([x0(i) y0(i)] - [xs ys]);
    % Distance between secondary sources and listener
    R = norm([x0(i) y0(i)] - [X Y]);
    % Time between secondary sources and virtual source
    t2(i) = R2/c;
    % === Time, in which pre-echos occur ===
    t(i) = R/c;
end

% Adjust the time, so the virtual source arrives at time 0 at the listener
% position and use only the minimum time (focused sources).
t = min(t-1.*t2 - norm([xs ys]-[X Y])/c);




