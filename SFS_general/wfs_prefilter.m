function hpre = wfs_prefilter(conf)
%WFS_PREFILTER creates a pre-equalization filter for WFS
%
%   Usage: hpre = wfs_prefilter([conf])
%
%   Input parameters:
%       conf - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       hpre - pre-equalization filter
%
%   WFS_PREFILTER(conf) calculates a sqrt(j k) pre-equalization filter for
%   Wave Field Synthesis (from conf.hpreflow to conf.hprefhigh, see SFS_config).
%
%   see also: wfs_preequalization, SFS_config, brs_wfs_25d

%*****************************************************************************
% Copyright (c) 2010-2012 Quality & Usability Lab                            *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
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
flow = conf.hpreflow;       % Lower frequency limit of preequalization
                            % filter (= frequency when subwoofer is active)
fhigh = conf.hprefhigh;     % Upper frequency limit of preequalization
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
% Apply sqrt(2*pi*f)/sqrt(2*pi*fhigh) filter for idxf < idxfhigh
H(1:idxfhigh) = sqrt(2*pi*f(1:idxfhigh))./sqrt(2*pi*fhigh);
% Set the response for idxf < idxflow to the value at idxflow
H(1:idxflow) = H(idxflow)*ones(1,idxflow);

% Compute filter
hpre = firls(Nfilt,2*f/fs,H);

% Truncate length to power of 2
hpre = hpre(1:end-1);
