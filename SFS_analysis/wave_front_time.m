function t = wave_front_time(X,xs,L,conf)
%WAVE_FRONT_TIME time of occurence of first echo for a linear WFS array
%   Usage: wave_front_time(X,xs,L,conf)
%          wave_front_time(X,xs,L)
%
%   Input parameters:
%       X       - listener position (m)
%       xs      - position of the virtual source (m)
%       L       - length of the linear loudspeaker array
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   WAVE_FRONT_TIME(X,xs,L) calculates the time of occuring of the first
%   echo (focused sources) resp. the last echo (virtual point source) at the
%   listener position X for a given point source location xs and a 
%   linear WFS loudspeaker array with a length of L.
%
%   see also: wave_front_direction, plot_wave_front_times
%

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input parameters ===================================
nargmin = 5;
nargmax = 6;
error(nargchk(nargmin,nargmax,nargin));
[X,xs] = position_vector(X,xs);
isargpositivescalar(L);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
% Loudspeaker distance
dx0 = conf.dx0;
% Array center position
X0 = conf.X0;
% Speed of sound
c = conf.c;
% Bandwidth of Dirac pulse
%f = fix(fs/2);
%k = 2*pi*f/c;


%% ===== Variables ======================================================
% Loudspeaker positions
x0 = secondary_source_positions(L,conf);
% Number of loudspeakers
nls = size(x0,1);


%% ===== Calculate a time axis ==========================================

% Geometry
%          x0(ii,:)                  X0
% x-axis <-^--^--^--^--^--^--^--^--^-|-^--^--^--^--^--^--^--^--^--
%             |                      |
%         R2 |  |                    |
%           |     | R                |      
%          x        |                |
%         xs          |              |
%                       O            |
%                       X            |
%                                    v
%                                  y-axis
%
% Distance between secondary source and listener: R
% Distance between secondary source and virtual source: R2

% Calculate arrival times at the virtual source position of the waves 
% emitted by the secondary sources
t2 = zeros(nls,1);
t = zeros(nls,1);
for i = 1:nls
    % Time between secondary sources and virtual source
    t2(ii) = norm(x0(ii,1:3)-xs)/c;
    % === Time, in which pre-echos occur ===
    t(ii) = norm(x0(ii,1:3)-X)/c;
end

% Adjust the time, so the virtual source arrives at time 0 at the listener
% position and use only the minimum time (focused sources).
t = min(t-1.*t2 - norm(xs-X)/c);
