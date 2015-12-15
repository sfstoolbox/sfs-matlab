function [hpreflow, hprefhigh] = localwfs_findhpref(X, phi, xs, src, conf)

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
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
nargmax = 5;
narginchk(nargmin,nargmax);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end

if conf.debug
    isargposition(X);
    isargxs(xs);
    isargscalar(phi);
    isargchar(src);
end

%% ===== Configuration ==================================================
conf.plot.useplot = false;  % disable plotting in easyfft
conf.ir.usehcomp = false;
conf.wfs.usehpre = false;     % no prefilter
conf.localsfs.wfs = conf.wfs;
%% ===== Variables ======================================================
N = conf.N;
irs = dummy_irs(N, conf);   % Impulse responses
fs = conf.fs;               % Sampling rate
dimension = conf.dimension; % dimensionality

%% ===== Computation ====================================================
% Compute impulse response/amplitude spectrum without prefilter
ir = ir_localwfs(X, phi, xs, src, irs, conf);
[H,~,f]=easyfft(ir(:,1),conf);

H = H./H(1);  % Normalize amplitude spectrum with H(f=0Hz)

% Model of local WFS spectrum without prefilter:
%   ^
% 1_| ______flow        
%   |       \          
%   |        \   
%   |         \_______
%   |          fhigh
%   -------------------------> f

if strcmp('2.5D',dimension)  
    % Expected slope: 6dB per frequency-doubling
    % Find 6dB cut-off frequency
    flowidx = find(H <= 1/2, 1, 'first');
    hpreflow = f(flowidx)/2;

elseif strcmp('3D',dimension) || strcmp('2D',dimension)
    % Expected slope: 12dB per frequency-doubling
    % Find 12dB cut-off frequency 
    flowidx = find(H <= 1/4, 1, 'first');
    hpreflow = f(flowidx)/4;  
    
else
    error('%s: %s is not a valid conf.dimension entry',upper(mfilename));
end

% approximated slope beginning at hpreflow
Hslope = hpreflow./f(flowidx:end);
% mean of H(f) evaluated from f to fs/2
Hmean = cumsum(H(end:-1:flowidx));  % cumulative sum
Hmean = Hmean./(1:length(Hmean)).';  % cumulative mean
Hmean = fliplr(Hmean);
% fhighidx is the frequency where both functions intersect the first time
fhighidx = flowidx - 1 + find( Hslope <= Hmean, 1, 'first');

if isempty(fhighidx)
  hprefhigh = fs/2;
else
  hprefhigh = f(fhighidx);
end 
