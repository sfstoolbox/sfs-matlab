function [alpha,a,t] = wave_front_direction(X,phi,xs,L,src,conf)
%WAVE_FRONT_DIRECTION returns direction, amplitude and time of the single wave
%   fronts for WFS
%
%   Usage: [alpha,a,t] = wave_front_direction(X,phi,xs,L,conf)
%          [alpha,a,t] = wave_front_direction(X,phi,xs,L)
%
%   Input parameters:
%       X,phi   - listener position and direction (m),(rad)
%       xs      - virtual source position (m)
%       L       - length of the linear loudspeaker array (m)
%       src     - used source model:
%                   'pw' - plane wave
%                   'ps' - point source
%                   'fs' - focused source
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       alpha   - angle of incident for every echo (rad) 
%       a       - amplitudes of the echos
%       t       - time of the wave fronts (s)
%
%   WAVE_FRONT_DIRECTION(X,phi,xs,src,L) calculates the direction of the single wave
%   fronts (due to aliasing artifacts) arriving from the loudspeakers for a
%   WFS array at the given listener position X for the given virtual source
%   xs and given array length L.
%
%   see also: brs_wfs_25d
%

% AUTHOR: Hagen Wierstorf, Sascha Spors
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input parameters ====================================
nargmin = 5;
nargmax = 6;
error(nargchk(nargmin,nargmax,nargin));
[X,xs] = position_vector(X,xs);
isargscalar(phi);
isargpositivescalar(L);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
% Loudspeaker distance
dx0 = conf.dx0;
% Speed of sound
c = conf.c;
% use plotting?
useplot = conf.useplot;
% Use gnuplot?
usegnuplot = conf.plot.usegnuplot;


%% ===== Variables =======================================================
phi = correct_azimuth(phi);
% name of the data file to store the result
outfile = sprintf('direction_L%i_xs%i_ys%i_X%.1f_Y%.1f.txt', ...
    L,xs(1),xs(2),X(1),X(2));
outfiledB = sprintf('direction_L%i_xs%i_ys%i_X%.1f_Y%.1f_dB.txt',...
                    L,xs(1),xs(2),X(1),X(2));

% Loudspeaker positions
x0 = secondary_source_positions(L,conf);
% Number of loudspeaker
nls = number_of_loudspeaker(L,conf);

% === Design tapering window ===
% See SFS_config if it is applied
win = tapwin(L,conf);

%% ===== Calculate direction of the echos ===============================

% Geometry
%           x0(ii,:)                 X0
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
% Distance between x0 and listener: R
% Distance between x0 and virtual source: R2

alpha = zeros(nls,1);
a = zeros(nls,1);
v = zeros(nls,2);
vdB = zeros(nls,2);
t = zeros(nls,1);
for ii = 1:nls

    % Get time delay and weighting factor for a single echo
    [weight,delay] = ...
        driving_function_imp_wfs_25d(x0(ii,:),xs,src,conf);

    % === Time, in which pre-echos occur ===
    t(ii) = norm(X-x0(ii,1:3))/c + delay - norm(X-xs)/c;

    % === Applying tapering window to the amplitude
    a(ii) = weight * win(ii);

    % === Direction of the wave fronts (in radian) ===
    % Angle between listener and secondary source (-pi < alpha <= pi, 
    % without phi)
    % Note: phi is the orientation of the listener (see first graph)
    alpha(ii) = cart2sph(x0(n,1)-X(1),x0(n,2)-X(2),0) - phi;

    % === Direction and amplitude in vector notation ===
    % Rotation matrix (see: http://en.wikipedia.org/wiki/Rotation_matrix)
    % RM = [cos(alpha(i)) -sin(alpha(i)); ...
    %       sin(alpha(i)) cos(alpha(i))];
    RM = rotation_matrix(alpha(ii),'counterclockwise');
    % Vector notation of angle and amplitude (in x,y coordinates)
    v(ii,:) = a(ii) .* (RM * [0 1]');
    %                 \____  ____/
    %                      \/
    %           unit vector in direction
    %           of the given loudspeaker
    % In dB notation
    vdB(ii,:) = (20*log10(a(ii))+100) .* (RM * [0 1]');
end


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
    alpha_ds = cart2sph(xs(1)-X(1),xs(2)-X(2),0) - phi;
    RM_ds = rotation_matrix(alpha_ds,'counterclockwise');
    X_ds = A .* (RM_ds * [0 1]');
    X_dsdB = (20*log10(A)+100) .* (RM_ds * [0 1]');

    % Position of the listener
    Xn = X(1).*ones(nls,1);
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
    fprintf(fid,'\n\n%f\t0\t%f\t%f\t%f',X(1),X_ds(1),X_ds(2),A);
    fprintf(fiddB,'\n\n%f\t0\t%f\t%f\t%f',X(1),X_dsdB(1),X_dsdB(2),20*log10(A)+100);
end
