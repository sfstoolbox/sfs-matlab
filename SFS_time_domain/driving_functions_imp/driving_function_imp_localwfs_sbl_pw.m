function [d, delay_offset] = driving_function_imp_localwfs_sbl_pw(x0, nk, conf)

%*****************************************************************************
% Copyright (c) 2016      Fiete Winter                                       *
%                         Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock, Germany  *
%                                                                            *
% This file is part of the supplementary material for Fiete Winter's         *
% scientific work and publications                                           *
%                                                                            *
% You can redistribute the material and/or modify it  under the terms of the *
% GNU  General  Public  License as published by the Free Software Foundation *
% , either version 3 of the License,  or (at your option) any later version. *
%                                                                            *
% This Material is distributed in the hope that it will be useful, but       *
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
% or FITNESS FOR A PARTICULAR PURPOSE.                                       *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy of the GNU General Public License along   *
% with this program. If not, see <http://www.gnu.org/licenses/>.             *
%                                                                            *
% http://github.com/fietew/publications           fiete.winter@uni-rostock.de*
%*****************************************************************************

%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 3;
narginchk(nargmin,nargmax);

%% ===== Configuration ========================================================
N0 = size(x0,1);
xref = conf.xref; 
fs = conf.fs;
% maximum order of circular basis expansion of sound field
if isempty(conf.localsfs.sbl.order)
    Nce = nfchoa_order(N0,conf);
else
    Nce = conf.localsfs.sbl.order;
end
% resolution of plane wave decomposition
if isempty(conf.localsfs.sbl.Npw)
  Npw = 2*ceil(2*pi*0.9*fs/conf.c*conf.secondary_sources.size/2);
else
  Npw = conf.localsfs.sbl.Npw;
end

wfsconf = conf;
wfsconf.wfs = conf.localsfs.wfs;
wfsconf.driving_functions = conf.localsfs.sbl.driving_functions;

%% ===== Computation ==========================================================
% circular expansion coefficients
[pm, delay_circexp] = circexp_imp_pw(nk, Nce, xref, conf);
% plane wave decomposition
ppwd = pwd_imp_circexp(pm, Npw);
% driving signal
[d, delay_lwfs] = driving_function_imp_wfs_pwd(x0, ppwd, xref, wfsconf);
% delay
delay_offset = delay_lwfs+delay_circexp;
