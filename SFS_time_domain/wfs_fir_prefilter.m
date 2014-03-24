function hpre = wfs_fir_prefilter(conf)
%WFS_FIR_PREFILTER creates a pre-equalization filter for WFS
%
%   Usage: hpre = wfs_fir_prefilter([conf])
%
%   Input parameters:
%       conf - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       hpre - pre-equalization filter
%
%   WFS_FIR_PREFILTER(conf) calculates a sqrt(j k) pre-equalization filter for
%   Wave Field Synthesis (from conf.wfs.hpreflow to conf.wfs.hprefhigh,
%   see SFS_config).
%
%   References:
%       H. Wierstorf (2014) - "Perceptual Assessment of Sound Field Synthesis",
%       PhD thesis, Tu Berlin
%
%   see also: wfs_preequalization, wfs_iir_prefilter, sound_field_imp_wfs, ir_wfs

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
nargmin = 0;
nargmax = 1;
narginchk(nargmin,nargmax);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
fs = conf.fs;               % Sampling rate
dimension = conf.dimension; % dimensionality
flow = conf.wfs.hpreflow;   % Lower frequency limit of preequalization
                            % filter (= frequency when subwoofer is active)
fhigh = conf.wfs.hprefhigh; % Upper frequency limit of preequalization
                            % filter (= aliasing frequency of system)


%% ===== Variables ======================================================

% Number of coefficients for filter
Nfilt=128;
% Frequency axis
f = linspace(0,fs/2,fs/10);
% Find indices for frequencies in f smaller and nearest to fhigh and flow
idxfhigh = max(find(f<fhigh));
idxflow = max(find(f<flow));
% Initialize response
H = ones(1,length(f));


%% ===== Computation ====================================================

% Desired response
%   ^
% 1_|          fhigh_______
%   |            /
%   |        /
%   | ___/
%   |  flow
%   -------------------------> f
%
% Pre-equilization filter from flow to fhigh
if strcmp('2.5D',dimension)
    %           _______
    %  H(f) = \|f/fhigh, for flow<=f<=fhigh
    %
    %  compare Wierstorf (2014), p. 25 (2.46)
    %
    H(idxflow:idxfhigh) = sqrt(f(idxflow:idxfhigh)./fhigh);
elseif strcmp('3D',dimension) || strcmp('2D',dimension)
    %         
    %  H(f) = f/fhigh, for flow<=f<=fhigh
    %
    %  compare Wierstorf (2014), p. 25 (2.45)
    %
    H(idxflow:idxfhigh) = f(idxflow:idxfhigh)./fhigh;
else
    error('%s: %s is not a valid conf.dimension entry',upper(mfilename));
end
% % Set the response for idxf < idxflow to the value at idxflow
H(1:idxflow) = H(idxflow)*ones(1,idxflow);

% Compute filter
hpre = firls(Nfilt,2*f/fs,H);

% Truncate length to power of 2
hpre = hpre(1:end-1);
