function gp_save_matrix(file,x,y,M)
% GP_SAVE_MATRIX save x,y,M in the matrix binary format of Gnuplot
%
%   Usage: gp_save_matrix(file,x,y,M)
%
%   Input parameters:
%       file    - filename of the data file
%       x       - x axis values (vector or matrix)
%       y       - y axis values (vector or matrix)
%       M       - matrix data size(M)=y,x
%
%   GP_SAVE_MATRIX(file,x,y,M) saves the values of x,y and M in a binary matrix
%   format useable by Gnuplot, see
%   http://www.gnuplotting.org/manpage-gnuplot-4-6/#x1-372000III for details.
%   If x,y are provided as matrices as, for example, they are returned by any
%   sound_field_*() function, x(1,:) and y(:,1) will be used.
%
%   See also: gp_save

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


%% ===== Checking of input  parameters ==================================
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargchar(file);
isargnumeric(x,y);
isargmatrix(M);


%% ===== Computation =====================================================
% Gnuplot can not handle arbritrary grids, so we extract the first row and first
% column instead
if is_dim_custom(x), x = x(1,:); end
if is_dim_custom(y), y = y(:,1); end
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
