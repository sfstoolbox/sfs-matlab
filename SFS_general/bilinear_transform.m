function [b,a] = bilinear_transform(sos,conf)
%BILINEAR_TRANSFORM returns the second-order section as filter coefficients
%
%   Usage: [b,a] = bilinear_transform(sos,[conf])
%
%   Input parameters:
%       sos     - second-order section representation
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       b,a     - filter coefficients as cells
%
%   BILINEAR_TRANSFORM(sos) transforms the second-order section representation
%   of a system to time discrete filter coefficients.
%
%   see also: driving_function_imp_nfchoa, zp2sos

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


%% ===== Checking input parameters =======================================
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
isargmatrix(sos);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
fs = conf.fs;


%% ===== Computation =====================================================
a = cell(size(sos,1),1);
b = a;
for n=1:size(sos,1)
    if isoctave
        [bz,az] = bilinear(sos(n,1:3),sos(n,4:6),1/fs);
    else
        [bz,az] = bilinear(sos(n,1:3),sos(n,4:6),fs);
    end 
    b{n} = bz;
    a{n} = az;
end
