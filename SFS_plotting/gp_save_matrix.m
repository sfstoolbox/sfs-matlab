function gp_save_matrix(file,x,y,M)
% GP_SAVE_MATRIX save x,y,M in the matrix binary format of Gnuplot
%
%   Usage: gp_save_matrix(file,x,y,M)
%
%   Input parameters:
%       file    - filename of the data file
%       x       - x axis values
%       y       - y axis values
%       M       - matrix data size(M)=y,x
%
%   GP_SAVE_MATRIX(file,x,y,M) saves the values of x,y and M in a binary matrix
%   format useable by Gnuplot, see
%   http://www.gnuplotting.org/manpage-gnuplot-4-6/#x1-372000III for details.
%
%   see also: gp_save

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


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargchar(file);
isargvector(x,y);
isargmatrix(M);


%% ===== Computation =====================================================

% Check if the data has the right format
[ly,lx] = size(M);
if lx~=length(x) || ly~=length(y)
    error('%s: size(M) has to be y,x!',upper(mfilename));
end

% Create matrix to store in the file
MS = zeros(length(x)+1,length(y)+1);
MS(1,1) = length(x);
MS(1,2:end) = y;
MS(2:end,1) = x;
MS(2:end,2:end) = M';

% Write data into the file
fid = fopen(file,'w');
fwrite(fid,MS,'float');
fclose(fid);
