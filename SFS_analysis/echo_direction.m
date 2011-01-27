function [alpha,a,t] = echo_direction(X,Y,phi,xs,ys,L,conf)
%ECHO_DIRECTION directions of the echos for a linear WFS array
%   Usage: [alpha,a,t] = echo_direction(X,Y,phi,xs,ys,L,conf)
%          [alpha,a,t] = echo_direction(X,Y,phi,xs,ys,L)
%
%   Input parameters:
%       X,Y,phi - listener position and direction
%       xs,ys   - virtual source position
%       L       - length of the linear loudspeaker array
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       alpha   - angle of incident for every echo (in °) 
%       a       - amplitudes of the echos
%       t       - time of the echos (s)
%
%   ECHO_DIRECTION(X,Y,phi,xs,ys,L) calculates the direction of the echos 
%   (due to aliasing artefacts) arriving from the loudspeakers for a linear
%   WFS array (length L) at the given listener position [X Y] for the given 
%   virtual source [xs ys] and given array length L.
%
%   see also: wfs_brs, wfs_amplitude_linear
%

% AUTHOR: Hagen Wierstorf, Sascha Spors

% TODO: 
%   * check if we have the right start time of the virtual source


%% ===== Checking of input  parameters ==================================

if nargchk(6,7,nargin)
    error(['Wrong number of args.',...
           'Usage: [alpha,a,t] = echo_direction(X,Y,phi,xs,ys,L,conf)']);
end

if ~isnumeric(X) || ~isscalar(X)
    error('%s: X has to be a scalar!',upper(mfilename));
end

if ~isnumeric(Y) || ~isscalar(Y)
    error('%s: Y has to be a scalar!',upper(mfilename));
end

if ~isnumeric(phi) || ~isscalar(phi)
    error('%s: phi has to be a scalar!',upper(mfilename));
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

if nargin<7
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

% use plotting?
useplot = conf.useplot;
% Use gnuplot?
usegnuplot = conf.usegnuplot;


%% ===== Variables ======================================================

% name of the data file to store the result
outfile = sprintf('direction_L%i_xs%i_ys%i_X%.1f_Y%.1f.txt',L,xs,ys,X,Y);
outfiledB = sprintf('direction_L%i_xs%i_ys%i_X%.1f_Y%.1f_dB.txt',...
                    L,xs,ys,X,Y);

% Number of loudspeaker (round towards plus infinity)
nLS = fix(L/LSdist)+1;

% Loudspeaker positions
[LSpos,LSdir] = LSpos_linear(X0,Y0,(nLS-1)*LSdist,nLS);

% === Design tapering window ===
% See SFS_config if it is applied
win = tapwin(L,conf)


%% ===== Calculate a time axis ==========================================

% Geometry
%           LSpos(:,i)            [X0 Y0]
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
% Distance between LSpos and listener: R
% Distance between LSpos and virtual source: R2

% Calculate arrival times at the virtual source position of the waves 
% emitted by the second sources
LStimes = zeros(nLS,1);
for i = 1:nLS
    % Distance between LSpos and virtual source
    R2 = norm([xs; ys] - LSpos(:,i));
    % Time between LSpos and virtual source
    LStimes(i) = R2/c;
end

% Start time of virtual source
% First loudspeaker signal for virtual source, last signal for focused
% sources.
if ys > 0
    xstime = max(LStimes);
else
    xstime = min(LStimes);
end

% Change order of the time vector, so the virtual source arrives at time 0 
% at the listener position.
t = -1.*LStimes - norm([xs ys]-[X Y])/c;


%% ===== Calculate direction of the echos ===============================

% Geometry
%           LSpos(:,i)            [X0,Y0]
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
% Distance between LSpos and listener: R
% Distance between LSpos and virtual source: R2

alpha = zeros(nLS,1);
a = zeros(nLS,1);
v = zeros(nLS,2);
vdB = zeros(nLS,2);
for i = 1:nLS
    
    % Distance between LSpos and listener (see Geometry)
    R = norm([X Y] - LSpos(:,i)');
    
    % === Time, in which pre-echos occur ===
    t(i) = R/c + t(i);
    
    % === Direction of the echos (in radian) ===
    % Vector from listener position to given loudspeaker position (cp. R)
    % NOTE: X - LSpos(1,n) gives a negative value for loudspeaker to the
    % left of the listener, therefor -dx1 is needed to get the right angle.
    dx = [X Y]' - LSpos(:,i);
    % Angle between listener and secondary source (-pi < alpha <= pi, 
    % without phi)
    % Note: phi is the orientation of the listener (see first graph)
    alpha(i) = atan2(-dx(1),dx(2)) - phi/180*pi;

    % === Amplitude factor ===
    % Use linear amplitude factor (see Spors et al. (2008)) and apply the
    % tapering window
    a(i) = wfs_amplitude_linear(LSpos(1,i),LSpos(2,i),X,Y,xs,ys) * win(i);

    % === Direction and amplitude in vector notation ===
    % Rotation matrix (see: http://en.wikipedia.org/wiki/Rotation_matrix)
    %RM = [cos(alpha(i)) -sin(alpha(i)); ...
    %      sin(alpha(i)) cos(alpha(i))];
    RM = rotation_matrix(alpha(i),'counterclockwise');
    % Vector notation of angle and amplitude (in x,y coordinates)
    v(i,:) = a(i) .* (RM * [0 1]');
    %                \____  ____/
    %                     \/
    %           unit vector in direction
    %           of the given loudspeaker
    % In dB notation
    vdB(i,:) = (20*log10(a(i))+100) .* (RM * [0 1]');
end

alpha = alpha./pi*180;

%% ===== Plotting =======================================================
if(useplot)
    % Plot the amplitude (dB) over time
    figure; plot(t,20*log10(a)+100);
end


%% ===== Generate data structures for plotting with gnuplot =============
if(usegnuplot)
    
    % === Amplitude and direction of the direct sound from the virtual
    % source ===
    % Calculate amplitude for the virtual source (which arrives per 
    % definition at t = 0
    A = sum(a);
    % Calculate the direction of the virtual source pulse (rad)
    alpha_ds = atan2(-(X-xs),Y-ys) - phi/180*pi;
    RM_ds = rotation_matrix(alpha_ds,'counterclockwise');
    X_ds = A .* (RM_ds * [0 1]');
    X_dsdB = (20*log10(A)+100) .* (RM_ds * [0 1]');
    
    % Position of the listener
    Xn = X.*ones(nLS,1);
    % Open file for storing results
    fid = fopen(outfile,'w');
    fiddB = fopen(outfiledB,'w');
    % Angle and amplitude in vector notation
    % Generate matrix for gnuplot plot
    gnuM = [Xn t v(:,1) v(:,2)];
    % Directions weighted by the amplitude of the echos
    %gnudB = [ Xn t sign(v(:,1)).*(20.*log10(abs(v(:,1)+eps))+100) sign(v(:,2)).*(20.*log10(abs(v(:,2)+eps))+100) ];
    gnudB = [Xn t vdB(:,1) vdB(:,2)];
    % Store data
    %save(outfile,'gnuM','-ascii');
    fprintf(fid,'#X\tt\tx\ty\n');
    fprintf(fiddB,'#X\tt\tx\ty\tA\n');
    fprintf(fid,'%f\t%f\t%f\t%f\n',gnuM');
    fprintf(fiddB,'%f\t%f\t%f\t%f\t%f\n',[gnudB'; (20*log10(a')+100)]);
    fprintf(fid,'\n\n%f\t0\t%f\t%f\t%f',X,X_ds(1),X_ds(2),A);
    fprintf(fiddB,'\n\n%f\t0\t%f\t%f\t%f',X,X_dsdB(1),X_dsdB(2),20*log10(A)+100);
end
