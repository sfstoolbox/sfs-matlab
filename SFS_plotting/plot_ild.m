function plot_ild(ild,phi)
%PLOT_ILD plots the given ILD values
%
%   Usage: plot_ild(ild,[phi])
%
%   Input options:
%       ild -   vector with given ILD values (e.g. created with
%               interaural_level_difference)
%       phi -   corresponding angles (optional, default: -180°..179°)
%
%   PLOT_ILD(ild,phi) creates a figure with the given ILD values and add
%   corresponding labels etc.
%
%   see also: interaural_level_difference, plot_itd

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

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 2;
error(nargchk(nargmin,nargmax,nargin));
isargvector(ild);
if nargin==nargmax
    isargvector(phi);
    if length(ild)~=length(phi)
        error('%s: phi has to have the same length as ild!',upper(mfilename));
    end
else
    phi = -180:1:179;
end


%% ===== Plotting ========================================================
figure;
plot(phi,ild);
axis([phi(1),phi(end),-30,30]);
xlabel('phi (°)');
ylabel('ILD (dB)')
