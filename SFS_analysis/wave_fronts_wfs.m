function varargout = wave_fronts_wfs(X,phi,xs,src,gnuplot,conf)
%WAVE_FRONTS_WFS returns direction, amplitude and time of the single wave
%   fronts for WFS
%
%   Usage: [alpha,a,t] = wave_fronts_wfs(X,phi,xs,src,[gnuplot],conf)
%
%   Input parameters:
%       X,phi   - listener position and direction / m, rad
%       xs      - virtual source position / m
%       src     - used source model:
%                   'pw' - plane wave
%                   'ps' - point source
%                   'fs' - focused source
%       gnuplot - boolean indicator if the data should be stored in to output
%                 files, called direction_nls*.txt for later plotting with
%                 gnuplot, default: false
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       alpha   - angle of incident for every echo / rad
%       a       - amplitudes of the echos
%       t       - time of the wave fronts / s
%
%   WAVE_FRONTS_WFS(X,phi,xs,src,conf) calculates the direction of the single
%   wave fronts (due to aliasing artifacts) arriving from the loudspeakers for a
%   WFS array at the given listener position X for the given virtual source
%   xs.
%
%   Refrences:
%   H. Wierstorf, A. Raake, M. Geier, S. Spors (2013) - Perception of Focused
%   Sources in Wave Field Synthesis, J. Audio Eng. Soc. 61.1, p. 5-16
%
%   See also: ir_wfs, driving_function_imp_wfs

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input parameters ====================================
nargmin = 5;
nargmax = 6;
narginchk(nargmin,nargmax);
isargposition(X);
isargxs(xs),
isargscalar(phi);
if nargin<nargmax
    conf = gnuplot;
    gnuplot = false;
end
isargstruct(conf);


%% ===== Configuration ===================================================
c = conf.c;
useplot = conf.plot.useplot;


%% ===== Variables =======================================================
phi = correct_azimuth(phi);
% Loudspeaker positions
x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,xs,src);
x0 = secondary_source_tapering(x0,conf);
% Number of loudspeaker
nls = size(x0,1);

%% ===== Calculate direction of the echos ===============================

% Geometry
%           x0(ii,:)                 X0
% x-axis <-^--^--^--^--^--^--^--^--^-|-^--^--^--^--^--^--^--^--^--
%             |                      |
%            |  |                    |
%           |     |                  |
%          x        |                |
%         xs          |              |
%                       O            |
%                       X            |
%                                    v
%                                  y-axis
%

% Get time delay and weighting factor for a single echo
[~,delay,a] = driving_function_imp_wfs(x0,xs,src,conf);
% Time, in which pre-echos occur:
% time_from_secondary_source_to_listener + delay
t = vector_norm(bsxfun(@minus,X,x0(:,1:3)),2)./c + delay -norm(X-xs(:,1:3))/c;
% Set start to t=0
t = t-min(t);
% Direction of the wave fronts / rad 
% angle between listener and secondary source (-pi < alpha <= pi,
% without phi)
% NOTE: phi is the orientation of the listener (see first graph)
[alpha_tmp,~,~] = cart2sph(x0(:,1)-X(1),x0(:,2)-X(2),x0(:,3)-X(3));
alpha = alpha_tmp - phi;

% Return values
if nargout>0, varargout{1}=alpha; end
if nargout>1, varargout{2}=a; end
if nargout>2, varargout{3}=t; end


%% ===== Plotting =======================================================
if (nargout==0 || useplot)
    % Plot the amplitude (dB) over time
    figure; plot(t,20*log10(a)+100,'xb');
    xlabel('t / s');
    ylabel('Amplitude / dB');
end


%% ===== Generate data structures for plotting with gnuplot =============
if gnuplot

    v = zeros(nls,3);
    vdB = zeros(nls,3);
    for ii = 1:nls
        % === Direction and amplitude in vector notation ===
        % Rotation matrix (see: http://en.wikipedia.org/wiki/Rotation_matrix)
        % RM = [cos(alpha) -sin(alpha); ...
        %       sin(alpha) cos(alpha)];
        RM = rotation_matrix(alpha(ii),3,'counterclockwise');
        % Vector notation of angle and amplitude (in x,y coordinates)
        v(ii,:) = a(ii) .* (RM * [0 1 0]');
        %                   \____  ____/
        %                        \/
        %             unit vector in direction
        %             of the given loudspeaker
        % In dB notation
        vdB(ii,:) = (20*log10(a(ii))+100) .* (RM * [0 1 0]');
    end

    % name of the data file to store the result
    outfile = sprintf('direction_nls%i_xs%.1f_ys%.1f_X%.1f_Y%.1f.txt', ...
        nls,xs(1),xs(2),X(1),X(2));
    outfiledB = sprintf('direction_nls%i_xs%.f1_ys%.1f_X%.1f_Y%.1f_dB.txt',...
                    nls,xs(1),xs(2),X(1),X(2));


    % == Amplitude and direction of the direct sound from the virtual source ==
    % Calculate amplitude for the virtual source (which arrives per
    % definition at t = 0
    A = sum(a);
    % Calculate the direction of the virtual source pulse (rad)
    [alpha_ds,~,~] = cart2sph(xs(1)-X(1),xs(2)-X(2),xs(3)-X(3));
    alpha_ds = alpha_ds - phi;
    RM_ds = rotation_matrix(alpha_ds,3,'counterclockwise');
    X_ds = A .* (RM_ds * [0 1 0]');
    X_dsdB = (20*log10(A)+100) .* (RM_ds * [0 1 0]');

    % Position of the listener
    Xn = X(1).*ones(nls,1);
    % Open file for storing results
    fid = fopen(outfile,'w');
    fiddB = fopen(outfiledB,'w');
    % Angle and amplitude in vector notation
    % Generate matrix for gnuplot plot
    gnuM = [Xn t v];
    % Directions weighted by the amplitude of the echos
    gnudB = [Xn t vdB];
    % Store data
    fprintf(fid,'#X\tt\tx\ty\tz\n');
    fprintf(fiddB,'#X\tt\tx\ty\tA\n');
    fprintf(fid,'%f\t%f\t%f\t%f\t%f\n',gnuM);
    fprintf(fiddB,'%f\t%f\t%f\t%f\t%f\t%f\n',[gnudB (20*log10(a)+100)]);
    fprintf(fid,'\n\n%f\t0\t%f\t%f\t%f',X(1),X_ds(1),X_ds(2),A);
    fprintf(fiddB,'\n\n%f\t0\t%f\t%f\t%f',X(1),X_dsdB(1),X_dsdB(2),20*log10(A)+100);
    fclose(fid);
    fclose(fiddB);
end
