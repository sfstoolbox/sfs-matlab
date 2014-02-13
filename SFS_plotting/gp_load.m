function [x,y,header] = gp_load(file)
%GP_LOAD load x,y from a text file
%
%   Usage: [x,y,header] = gp_load(file)
%
%   Input parameter:
%       file    - filename of the data file
%
%   Output paramter:
%       x       - x axis values
%       y       - y axis values [vector or matrix]
%       header  - header comment
%
%   GP_LOAD(file) load the values of x and y from a text, that was saved
%   in a Gnuplot compatible format

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

% FIXME: is this function doing what it should do?


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargfile(file);


%% ===== Computation =====================================================

% FIXME: use fopen etc to read. Handle header line!
% Read the data
if isoctave
    to_be_implemented;
else
    content = textread(file,'','delimiter','\t','commentstyle','shell');
end

x = content(:,1);
y = content(:,2:end);
header = '';
