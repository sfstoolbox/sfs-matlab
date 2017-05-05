function [d,delay_offset] = driving_function_imp_wfs_pwd(x0,ppwd,xq,conf)
%DRIVING_FUNCTION_IMP_WFS_PWD returns the WFS driving signal for a plane wave 
%expansion
%
%   Usage: [d,delay_offset] = driving_function_imp_wfs_pwd(x0,ppwd,[xq],conf)
%
%   Input parameters:
%       x0          - position and direction of the secondary source / m [N0x6]
%       ppwd        - plane wave coefficients [N x Npw]
%       xq          - centre of plane wave expansion, default = [0,0,0];
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       d             - driving function signals [N x N0]
%       delay_offset  - additional added delay, so you can correct it
%
%   See also: driving_function_imp_wfs_vss, driving_function_imp_localwfs_sbl

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Developers                             *
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
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargsecondarysource(x0);
isargmatrix(ppwd);
if nargin == nargmin
    conf = xq;
    xq = [0, 0, 0];
else
    isargxs(xq);
end
isargstruct(conf);


%% ===== Computation ====================================================
% create distribution of plane waves as virtual secondary sources
Npw = size(ppwd, 2);
phipw = (0:Npw-1).'*2*pi/Npw;
xv = [cos(phipw), sin(phipw)];  % [Npw x 2]
xv(:,3) = 0;
xv(:,4:6) = xv(:,1:3);
xv(:,7) = 1./Npw;  % apply integrations weights to pwd

% shift coordinates to expansion center
x0(:,1:3) = bsxfun(@minus,x0(:,1:3),xq);
conf.xref = [0 0 0];

[d,delay_offset] = driving_function_imp_wfs_vss(x0,xv,'pw',ppwd,conf);
