function [ir,x0] = ir_nfchoa(X,phi,xs,src,irs,conf)
%IR_NFCHOA Generate a IR for NFCHOA
%
%   Usage: ir = ir_nfchoa(X,phi,xs,src,irs,[conf])
%
%   Input parameters:
%       X       - listener position / m
%       phi     - listener direction [head orientation] / rad
%                 0 means the head is oriented towards the x-axis.
%       xs      - virtual source position [ys > Y0 => focused source] / m
%       src     - source type: 'pw' -plane wave
%                              'ps' - point source
%       irs     - IR data set for the secondary sources
%       conf    - optional configuration struct (see SFS_config) 
%
%   Output parameters:
%       ir      - impulse response for the desired HOA synthesis (nx2 matrix)
%       x0      - secondary sources
%
%   IR_NFCHOA(X,phi,xs,src,L,irs,conf) calculates a binaural room impulse
%   response for a virtual source at xs for a virtual NFCHOA array and a
%   listener located at X.
%
%   see also: ssr_brs_nfchoa, ir_nfchoa, ir_point_source, auralize_ir

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
nargmin = 5;
nargmax = 6;
narginchk(nargmin,nargmax);
if nargin==nargmax-1
    conf = SFS_config;
end
if conf.debug
    isargposition(X);
    isargxs(xs);
    isargscalar(phi);
    isargpositivescalar(L);
    isargchar(src);
    check_irs(irs);
end


%% ===== Variables ======================================================
% Loudspeaker positions
x0 = secondary_source_positions(conf);


%% ===== BRIR ===========================================================
% calculate driving function
d = driving_function_imp_nfchoa(x0,xs,src,conf);

% generate the impulse response for NFCHOA
ir = ir_generic(X,phi,x0,d,irs,conf);
