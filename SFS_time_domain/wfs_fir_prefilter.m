function hpre = wfs_fir_prefilter(conf)
%WFS_FIR_PREFILTER creates a pre-equalization filter for WFS
%
%   Usage: hpre = wfs_fir_prefilter(conf)
%
%   Input parameters:
%       conf - configuration struct (see SFS_config)
%
%   Output parameters:
%       hpre - pre-equalization filter
%
%   WFS_FIR_PREFILTER(conf) calculates a sqrt(j k) pre-equalization filter for
%   Wave Field Synthesis (from conf.wfs.hpreflow to conf.wfs.hprefhigh,
%   see SFS_config).
%
%   See also: wfs_preequalization, wfs_iir_prefilter, sound_field_imp_wfs, ir_wfs

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
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargstruct(conf);


%% ===== Configuration ==================================================
fs = conf.fs;                  % Sampling rate
c = conf.c;                    % Speed of sound
dimension = conf.dimension;    % Dimensionality
flow = conf.wfs.hpreflow;      % Lower frequency limit of preequalization
                               % filter (= frequency when subwoofer is active)
fhigh = conf.wfs.hprefhigh;    % Upper frequency limit of preequalization
                               % filter (= aliasing frequency of system)
Nfilt = conf.wfs.hpreFIRorder; % Number of coefficients for filter
if isodd(Nfilt)
    error(['%s: conf.wfs.hpreFIRorder == %i is not a valid filter order. ', ...
        'Must be an even integer.'],upper(mfilename),Nfilt);
end
%% ===== Variables ======================================================
% Frequency axis
f = linspace(0,fs/2,fs/10);
% Find indices for frequencies in f smaller and nearest to fhigh and flow
idxfhigh = max(find(f<fhigh));
idxflow = max(find(f<flow));
% Initialize response
H = ones(1,length(f));


%% ===== Computation ====================================================
% Desired response
%           ^
% H(fhigh) _|          fhigh_______
%           |            /
%           |        /
%  H(flow) _| ___/
%           |  flow
%           -------------------------> f
%
% Pre-equalization filter from flow to fhigh
if strcmp('2.5D',dimension)
    %            ______
    %           |2*pi*f
    %  H(f) = \ |------ for flow<=f<=fhigh
    %          \|  c
    %
    %  See http://sfstoolbox.org/#equation-h.wfs.2.5D
    %
    H(idxflow:idxfhigh) = sqrt(2*pi*f(idxflow:idxfhigh)/c);
    H(idxfhigh:end) = H(idxfhigh);
elseif strcmp('3D',dimension) || strcmp('2D',dimension)
    %
    %         2*pi*f
    %  H(f) = ------ for flow<=f<=fhigh
    %           c
    %
    %  See http://sfstoolbox.org/#equation-h.wfs
    %
    H(idxflow:idxfhigh) = 2*pi*f(idxflow:idxfhigh)/c;
    H(idxfhigh:end) = H(idxfhigh);
else
    error('%s: %s is not a valid conf.dimension entry',upper(mfilename),dimension);
end
% Set the response for idxf < idxflow to the value at idxflow
H(1:idxflow) = H(idxflow);

% Compute filter
hpre = firls(Nfilt,2*f/fs,H);

