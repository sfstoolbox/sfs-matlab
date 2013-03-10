function [x,y,P] = wave_field_mono_2d(X,Y,x0,D,f,conf)
%WAVE_FIELD_MONO_2D simulates a monofrequent wave field for line sources as
%secondary sources
%
%   Usage: [x,y,P] = wave_field_mono_2d(X,Y,x0,D,f,[conf])
%
%   Input parameters:
%       X           - [xmin,xmax]
%       Y           - [ymin,ymax]
%       x0          - secondary sources [n x 6]
%       D           - driving signals for the secondary sources [m x n]
%       f           - monochromatic frequency (Hz)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       x           - corresponding x axis
%       y           - corresponding y axis
%       P           - Simulated wave field
%
%   WAVE_FIELD_MONO(X,Y,x0,D,f,conf) simulates a wave field for the given
%   secondary sources, driven by the corresponding driving signals. Line
%   sources are applied as source models for the secondary sources. The
%   simulation is done for one frequency in the frequency domain, by
%   calculating the integral for P with a summation.
%   
%   To plot the result use plot_wavefield(x,y,P).
%
%   References:
%       
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: plot_wavefield, wave_field_mono_wfs_25d

%*****************************************************************************
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
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


%% ===== Checking of input  parameters ==================================
nargmin = 5;
nargmax = 6;
narginchk(nargmin,nargmax);
isargvector(X,Y,D);
isargsecondarysource(x0);
isargpositivescalar(f);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
if size(x0,1)~=length(D)
    error(['%s: The number of secondary sources (%i) and driving ', ...
        'signals (%i) does not correspond.'], ...
        upper(mfilename),size(x0,1),length(D));
end


%% ===== Configuration ==================================================
% Plotting result
useplot = conf.useplot;


%% ===== Computation ====================================================
% Create a x-y-grid
[xx,yy,x,y] = xy_grid(X,Y,conf);
% Initialize empty wave field
P = zeros(length(y),length(x));
% Integration over secondary source positions
for ii = 1:size(x0,1)

    % ====================================================================
    % Secondary source model G(x-x0,omega)
    % This is the model for the loudspeakers we apply. We use closed cabinet
    % loudspeakers and therefore point sources.
    G = line_source(xx,yy,x0(ii,1:3),f);

    % ====================================================================
    % Integration
    %              /
    % P(x,omega) = | D(x0,omega) G(x-x0,omega) dx0
    %              /
    %
    % see: Spors2009, Williams1993 p. 36
    P = P + D(ii).*G;

end

% === Scale signal (at xref) ===
P = norm_wave_field(P,x,y,conf);

% ===== Plotting =========================================================
if(useplot)
    plot_wavefield(x,y,P,x0,conf);
end
