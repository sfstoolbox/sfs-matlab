function hpre = wfs_iir_prefilter(conf)
%WFS_IIR_PREFILTER creates a minimum-phase IIR pre-equalization filter for WFS
%
%   Usage: hpre = wfs_iir_prefilter(conf)
%
%   Input parameters:
%       conf    - configuration struct (see SFS_config)
%
%       must at least include:
%           conf.fs = 44100;           - sampling frequency / Hz
%           conf.wfs.hpreflow = 200;   - lower shelving frequency for coupling
%                                        subwoofers and adapt the low frequency
%                                        to different array lengths / Hz
%           conf.wfs.hprefhigh = 1500; - higher shelving frequency to adapt to
%                                        actual aliasing frequeny / Hz
%           conf.wfs.hpreBandwidth_in_Oct = 2;  - bandwidth for the Lagrange
%                                                 interpolation region / octaves
%           conf.wfs.hpreIIRorder = 4; - desired IIR filter order
%
%   Output parameters:
%        hpre   - iir pre-equalization filter as a struct with the following
%                 entries:
%                   z,p,k       includes poles, zeros and gains equivalent to
%                   sos,g       includes a second-order section representation
%                               equivalent to sos2tf
%                   b,a         includes the filter coefficients (not suitable
%                               for higher order IIRs)
%                   max_dev_dB  maximum deviation in dB from desired shelving
%                               filter
%
%   Required Functions:
%   get_shelve_lagrange.m (included in SFS-toolbox)
%   fdesign (included in the Matlab Signal Processing Toolbox, requiring DSP System Design Toolbox)
%
%
%   WFS_IIR_PREFILTER(conf) calculates a sqrt(j k) pre-equalization filter
%   with high shelving characterstics for Wave Field Synthesis.
%   Note, this function does not work in Octave, use conf.wfs.hpretype='FIR'
%   instead.
%
%   For details see [Sch13]:
%   Frank Schultz, Vera Erbes, Sascha Spors, Stefan Weinzierl (2013):
%   "Derivation of IIR prefilters for soundfield synthesis using linear
%   secondary source distributions", In: Proc. of the
%   International Conference on Acoustics AIA-DAGA 2013, Merano, Italy,
%   18 - 21 March 2013, pages 2372-2375
%
%   See also: wfs_preequalization, wfs_fir_prefilter, driving_function_imp_wfs,
%       wfs_ir


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
% Revision: 07/02/2013 frank.schultz@uni-rostock.de initial development      *
%*****************************************************************************


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargstruct(conf);
% This function is not working in Octave at the moment.
if isoctave
    error(['%s: Not available under Octave, please use ', ...
           'conf.wfs.hpretype="FIR"'],upper(mfilename));
else
    v = ver;
    if ~any(strcmp('DSP System Toolbox', {v.Name}))
        error(['%s: The DSP System Toolbox is required for IIR filter, ', ...
               'please use conf.wfs.hpretype="FIR"'],upper(mfilename));
    end
end


%% ===== Configuration ==================================================
fs = conf.fs;               % Sampling rate
fsub = conf.wfs.hpreflow;       % Lower frequency limit of preequalization
                            % filter (= frequency when subwoofer is active)
falias = conf.wfs.hprefhigh;    % Upper frequency limit of preequalization
                            % filter (= aliasing frequency of system)

% Bandwidth in octaves for lagrange interpolation region
%at the moment only 0.5, 1,2,3 or 4
Bandwidth_in_Oct = conf.wfs.hpreBandwidth_in_Oct;
IIRorder = conf.wfs.hpreIIRorder;
debug = conf.debug;


%% ===== Variables ======================================================
NFFT = fs;                  %FFT length, note that we assume that fs and thus NFFT is EVEN!!!
df = fs/NFFT;               %FFT resolution is 1Hz!
f = (0:df:fs-df)';          %frequency vector


%% ===== Computation ====================================================
% Get +3dB/oct. ideal slope:
H_Pre3dB = 10*log10(f);
% Lagrange interpolation towards shelving filter:
H_Pre3dB_Shv_Lagrange = (get_shelve_lagrange(f,10.^(H_Pre3dB/20),1,fsub,1,falias,Bandwidth_in_Oct));
% Normalize to maximum, i.e. the flat amplitude response at high frequencies
% Note that for different falias the filter is louder or softer, which may be
% undesired, however we handle this for consistency with wfs_prefilter.m
% consider to normalize to 500 Hz instead in future:
H_Pre3dB_Shv_Lagrange = H_Pre3dB_Shv_Lagrange/max(abs(H_Pre3dB_Shv_Lagrange));
% Apply phase +pi/4, note that we actually don't need that when calling iirlpnorm.m
H_Pre3dB_Shv_Lagrange = H_Pre3dB_Shv_Lagrange*exp(+1i*pi/4);
% Prepare for iirlpnorm:
H = transpose(abs(H_Pre3dB_Shv_Lagrange(1:NFFT/2+1)));
F = transpose(f(1:NFFT/2+1)/(fs/2));
%****************
%IIR with LP-Norm
%****************
Nb = IIRorder; Na = Nb; %filter order for num/den
d = fdesign.arbmag('Nb,Na,F,A',Nb,Na,F,abs(H)); %we only consider abs, due to the desired minphase design
%help(d,'iirlpnorm')
% We call iirlpnorm with Matlab defaults, which have been checked regarding the
% consistency for versions 2010-2013
Hd_lpnorm = design(d, 'iirlpnorm');
if debug
isstable(Hd_lpnorm)
isminphase(Hd_lpnorm)
end
[hpre.z,hpre.p,hpre.k] = sos2zp(Hd_lpnorm.sosMatrix,Hd_lpnorm.ScaleValues);
% Don't use b,a for higher order IIR filters directly for stability reasons!
[hpre.b,hpre.a] =  sos2zp(Hd_lpnorm.sosMatrix,Hd_lpnorm.ScaleValues);
hpre.sos = Hd_lpnorm.sosMatrix;
hpre.g = Hd_lpnorm.ScaleValues;
% Check deviation in dB:
[H_Pre3dB_Shv_Lagrange_IIR] = freqz(Hd_lpnorm,f,fs);
% Plot:
if debug
    figure
    semilogx(f,20*log10(abs(H_Pre3dB_Shv_Lagrange)),'k','LineWidth',3), hold on
    semilogx(f,20*log10(abs(H_Pre3dB_Shv_Lagrange_IIR)),'r','LineWidth',2), hold off
    xlabel('f / Hz'), ylabel('A / dB'), title('IIR prefilter')
    axis([10 10000 -15 +3]), set(gca,'YTick',-15:3:+3)
    grid on
    legend('desired','minphase IIR realization')
end
