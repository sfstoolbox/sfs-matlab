function figsize(x,y,unit)
%FIGSIZE changes the size of a figure
%
%   Usage: figsize(x,y,unit)
%
%   Input options:
%       x,y         - x,y size of the figure
%       unit        - unit in which the size is given, can be one of the
%                     following: 'cm', 'px'
%
%   FIGSIZE(x,y,unit) sets the size of the last figure to x,y in the given unit.
%
%   see also: plot_wavefield

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
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
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input parameter =====================================
nargmin = 0;
nargmax = 3;
narginchk(nargmin,nargmax);
if nargin==0
    x = 10;
    y = 10;
    unit = 'cm';
end
isargpositivescalar(x,y)
isargchar(unit)


%% ===== Main ============================================================
% convert to centimeters
if strcmp('px',unit) || strcmp('pixel',unit)
    x = px2cm(x);
    y = px2cm(y);
elseif strcmp('inches',unit)
    x = in2cm(x);
    y = in2cm(y);
end

% Adjust the font size to match the figure size
%set_font_size(12);

% Set the figure dimensions
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[y,x]);
set(gcf,'PaperPosition',[0,0,x,y]);
