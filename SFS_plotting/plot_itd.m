function plot_itd(itd,phi)
%PLOT_ITD plots the given ITD values
%
%   Usage: plot_itd(itd,[phi])
%
%   Input options:
%       itd -   vector with given ITD values (e.g. crfeated with
%               interaural_time_difference)
%       phi -   corresponding angles (optional, default: -180°..179°)
%
%   PLOT_ITD(itd,phi) creates a figure with the given ITD values and add
%   corresponding labels etc.
%
%   see also: interaural_time_difference, plot_ild

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


% FIXME: is this function needed anymore?


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargvector(itd);
if nargin==nargmax
    isargvector(phi);
    if length(itd)~=length(phi)
        error('%s: phi has to have the same length as itd!',upper(mfilename));
    end
else
    phi = -180:1:179;
end


%% ===== Plotting ========================================================
figure;
plot(phi,1000*itd);
axis([phi(1),phi(end),-1,1]);
xlabel('phi (°)');
ylabel('ITD (ms)')
