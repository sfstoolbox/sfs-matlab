function [P,x0,win,weights] = wfs_3d(x,y,z,xs,src,f,L,conf)
%WFS_3D returns the sound pressure for 3D WFS at x,y,z
%
%   Usage: [P,x0,win] = wfs_3d(x,y,z,xs,src,f,L,[conf])
%
%   Input parameters:
%       x           - x position(s)
%       y           - y position(s)
%       z           - z position(s)
%       xs          - position of point source (m)
%       src         - source type of the virtual source
%                         'pw' - plane wave (xs is the direction of the
%                                plane wave in this case)
%                         'ps' - point source
%                         'fs' - focused source
%       f           - monochromatic frequency (Hz)
%       L           - array length (m)
%       conf        - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       P           - Simulated wave field
%       win         - tapering window (activity of loudspeaker)
%
%   WAVE_FIELD_MONO_WFS_3D(x,y,z,xs,L,f,src,conf) returns the sound pressure at
%   the point(s) (x,y,z) for the given source type (src) using a WFS 3
%   dimensional driving function in the temporal domain. This means by
%   calculating the integral for P with a summation.
%
%   References:
%       Spors2009 - Physical and Perceptual Properties of Focused Sources in
%           Wave Field Synthesis (AES127)
%       Williams1999 - Fourier Acoustics (Academic Press)
%
%   see also: wave_field_mono_wfs_25d

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


%% ===== Checking of input  parameters ==================================
nargmin = 7;
nargmax = 8;
error(nargchk(nargmin,nargmax,nargin));
isargmatrix(x,y);
xs = position_vector(xs);
isargpositivescalar(L,f);
isargchar(src);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ==
xref = conf.xref;


%% ===== Computation ====================================================
% Calculate the wave field in time-frequency domain
%
% Get the position of the loudspeakers and its activity
x0 = secondary_source_positions(L,conf);
x0 = secondary_source_selection(x0,xs,src,xref);
% Generate tapering window
win = tapering_window(x0(:,1:6),conf);
% Initialize empty wave field
% FIXME: it could be that length is not enough here and we need size(...)
P = zeros(length(y),length(x));
% Integration over secondary source positions

phi = zeros(1,size(x0,1));

if conf.debug
    D_plot =  zeros(1,size(x0,1));
end

for ii = 1:size(x0,1)

    % ====================================================================
    % Secondary source model G(x-x0,omega)
    % This is the model for the loudspeakers we apply. We use closed cabinet
    % loudspeakers and therefore point sources.
    G = point_source(x,y,z,x0(ii,1:3),f);

    % ====================================================================
    % Driving function D(x0,omega)
    D = driving_function_mono_wfs_3d(x0(ii,:),xs,src,f,conf);
    % ====================================================================
    % Integration
    %              /
    % P(x,omega) = | D(x0,omega) G(x-x0,omega) dx0
    %              /
    %
    % see: Spors2009, Williams1993 p. 36
    %
    % NOTE: win(ii) is the factor of the tapering window in order to have fewer
    % truncation artifacts. If you don't use a tapering window win(ii) will
    % always be one.
    
    P = P + win(ii)*D.*G;

    if conf.debug
        phi(1,ii) = atan2(x0(ii,2),x0(ii,1));
%           if phi(1,ii)<-2.9 || phi(1,ii)>-0.2
%               D_test = D*0.2;
%           else
%               D_test = D;
%           end
        D_plot(1,ii) = D;
%         D_plot2(1,ii) = D_test;
    end
    

%     P = P + win(ii)*D_test.*G;
    
end

    if conf.debug
%     save('DrivingFunctionWithoutWeights.mat','D_plot')  
    figure
    plot(phi(1,:),20*log10(abs(D_plot(1,:))),'k.');
    hold on 
%     plot(phi(1,:),20*log10(abs(D_plot2(1,:))),'r.');
    grid on
    xlabel('$\varphi$ / rad','interpreter','latex')
    ylabel('amplitude / dB')
    title('amplitude for all N loudspeakers')
    end