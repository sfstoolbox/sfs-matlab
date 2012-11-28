function [alpha,a,t] = wave_front_direction(X,phi,xs,src,L,conf)
%WAVE_FRONT_DIRECTION returns direction, amplitude and time of the single wave
%   fronts for WFS
%
%   Usage: [alpha,a,t] = wave_front_direction(X,phi,xs,src,L,[conf])
%
%   Input parameters:
%       X,phi   - listener position and direction (m),(rad)
%       xs      - virtual source position (m)
%       src     - used source model:
%                   'pw' - plane wave
%                   'ps' - point source
%                   'fs' - focused source
%       L       - length of the linear loudspeaker array (m)
%       conf    - optional configuration struct (see SFS_config)
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
%   see also: ir_wfs_25d

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************

% FIXME: see if this function still works correctly


%% ===== Checking of input parameters ====================================
nargmin = 5;
nargmax = 6;
narginchk(nargmin,nargmax);
[X,xs] = position_vector(X,xs);
isargscalar(phi);
isargpositivescalar(L);

if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
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
nls = size(x0,1);

% === Design tapering window ===
% See SFS_config if it is applied
win = tapering_window(x0,conf);

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
    [alpha_tmp,theta_tmp,r_tmp] = cart2sph(x0(ii,1)-X(1),x0(ii,2)-X(2),0);
    alpha(ii) = alpha_tmp - phi;

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
    [alpha_ds,theta_tmp,r_tmp] = cart2sph(xs(1)-X(1),xs(2)-X(2),0) - phi;
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
