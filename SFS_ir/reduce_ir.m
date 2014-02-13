function ir = reduce_ir(ir,fs,nsamples,conf)
%REDUCE_IR resamples and shortens a IR
%
%   Usage: ir = reduce_ir(ir,fs,nsamples,[conf])
%
%   Input parameters:
%       ir          - two channel IR signal
%       fs          - sampling rate of the target IR / Hz
%       nsamples    - length of the target IR
%       conf        - optional configuration struct (see SFS_config)
%
%   Output paramteres:
%       ir          - two channel IR signal
%
%   REDUCE_IR(ir,fs,nsamples,conf) shortens and resamples a given IR.
%   This can be useful for mobile phones.
%
%   see also: SFS_config, read_irs, intpol_ir, shorten_ir

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
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargpositivescalar(fs,nsamples);
if ~isnumeric(ir) || size(ir,2)~=2
    error('%s: ir has to be an IR with samples x 2 size.',upper(mfilename));
end
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
ofs = conf.fs;  % original fs
useplot = conf.plot.useplot;


%% ===== Computation ====================================================
% Resample HRIR
if ofs~=fs
    resamp_ir(:,1) = resample(ir(:,1),fs,ofs);
    resamp_ir(:,2) = resample(ir(:,2),fs,ofs);
else
    resamp_ir = ir;
end

% Window HRIR
win = hann_window(ceil(0.15*nsamples),ceil(0.10*nsamples),nsamples).^2;

% Find maximum of resampled HRIR
% Find maximum in each channel and calculate the mean of the index
[~,idx1] = max(abs(resamp_ir(:,1)));
[~,idx2] = max(abs(resamp_ir(:,2)));
idx = round((idx1+idx2)/2);

% Cut the HRIR around the maximum
% Leading zeros before idx
offset = 24;
short_ir(1:nsamples,1) = ...
    resamp_ir(idx-offset:idx+nsamples-offset-1,1) .* win;
short_ir(1:nsamples,2) = ...
    resamp_ir(idx-offset:idx+nsamples-offset-1,2) .* win;


%% ===== Plotting =======================================================

if(useplot)
    figure
    plot(resamp_ir(:,1),'-b'); hold on;
    plot(resamp_ir(:,2),'-r');
    figure
    plot(short_ir(:,1),'-b'); hold on;
    plot(short_ir(:,2),'r-');
end

