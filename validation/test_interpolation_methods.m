function status = test_interpolation_methods(modus)
%TEST_INTERPOLATION_METHODS demonstrates the difference between interpolation
%methods
%
%   Usage: status = test_interpolation_methods(modus)
%
%   Input parameters:
%       modus   - 0: numerical (not available)
%                 1: visual
%
%   Output parameters:
%       status  - true or false
%
%   TEST_INTERPOLATION_METHODS(modus) demonstrates the different interpolation
%   methods in the interpolate_ir function for shifted Dirac impulses and for
%   HRIRs of different directions

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


status = false;


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
if ~modus
    warning('%s: numerical modus not available.',upper(mfilename));
    status = true;
    return;
end


%% ===== Settings ========================================================
conf = struct;
conf.ir.useinterpolation = true;


%% ===== Main ============================================================
% Check if HRTF data set with low frequency correction is available,
% download otherwise
basepath = get_sfs_path();
hrtf = 'KEMAR_HRTFs_lfcorr.sofa';
hrtf_file = [basepath '/data/HRTFs/' hrtf];
if ~exist(hrtf_file,'file')
    disp('Download');
    url = ['https://github.com/spatialaudio/lf-corrected-kemar-hrtfs/', ...
           'blob/master/' hrtf '?raw=true'];
    download_file(url,hrtf_file);
end
% Load HRTF data set
hrtf = SOFAload(hrtf_file);


%% ===== Interpolate shifted Dirac impulses ==============================
% interpolation points
x0 = [1 3; 0 0; 0 0]/180*pi;
% weights (for target point: xs = [2; 0; 0]/180*pi )
weights = [.5 .5];
N = 50; %length impulse responses

% 1. Interpolate between Dirac impulses with one sample in between
h1 = zeros(1,N);
h1(3) = 1;
h2 = zeros(1,N);
h2(5) = 1;

irs1(1,1,:) = h1;
irs1(2,1,:) = h2;
irs1(:,2,:) = irs1(:,1,:); %redundant second channel

conf.ir.interpolationmethod = 'simple';
h_int1_simple = interpolate_ir(irs1,weights,conf);
conf.ir.interpolationmethod = 'freqdomain';
h_int1_fd = interpolate_ir(irs1,weights,conf);

% 2. Interpolate between neighbouring Dirac impulses at impulse reponse start
h1 = zeros(1,N);
h1(3) = 1;
h3 = zeros(1,N);
h3(4) = 1;

irs2(1,1,:) = h1;
irs2(2,1,:) = h3;
irs2(:,2,:) = irs2(:,1,:); %redundant second channel

conf.ir.interpolationmethod = 'simple';
h_int2_simple = interpolate_ir(irs2,weights,conf);
conf.ir.interpolationmethod = 'freqdomain';
h_int2_fd = interpolate_ir(irs2,weights,conf);

% 3. Interpolate between neighbouring Dirac impulses in middle of impulse response
h4 = zeros(1,N);
h4(25) = 1;
h5 = zeros(1,N);
h5(26) = 1;

irs3(1,1,:) = h4;
irs3(2,1,:) = h5;
irs3(:,2,:) = irs3(:,1,:);

conf.ir.interpolationmethod = 'simple';
h_int3_simple = interpolate_ir(irs3,weights,conf);
conf.ir.interpolationmethod = 'freqdomain';
h_int3_fd = interpolate_ir(irs3,weights,conf);

% Plots
% impulse responses
figure
    plot(0:N-1,squeeze(irs1(1,1,:)),'k'), hold on
    plot(0:N-1,squeeze(irs1(2,1,:)),'b')
    plot(0:N-1,squeeze(h_int1_simple(1,1,:)),'r')
    plot(0:N-1,squeeze(h_int1_fd(1,1,:)),'m')
    grid
    xlabel('samples'), ylabel('amplitude')
    legend('h_1','h_2','simple interp','freqdomain interp')
    title('Interpolation between Dirac impulses with one sample in between')

figure
    plot(0:N-1,squeeze(irs2(1,1,:)),'k'), hold on
    plot(0:N-1,squeeze(irs2(2,1,:)),'b')
    plot(0:N-1,squeeze(h_int2_simple(1,1,:)),'r')
    plot(0:N-1,squeeze(h_int2_fd(1,1,:)),'m')
    grid
    xlabel('samples'), ylabel('amplitude')
    legend('h_1','h_2','simple interp','freqdomain interp')
    title('Interpolation of neighbouring Dirac impulses at impulse response start')

figure
    plot(0:N-1,squeeze(irs3(1,1,:)),'k'), hold on
    plot(0:N-1,squeeze(irs3(2,1,:)),'b')
    plot(0:N-1,squeeze(h_int3_simple(1,1,:)),'r')
    plot(0:N-1,squeeze(h_int3_fd(1,1,:)),'m')
    grid
    xlabel('samples'), ylabel('amplitude')
    legend('h_1','h_2','simple interp','freqdomain interp')
    title('Interpolation of neighbouring Dirac impulses in middle of impulse response')


%% ===== Interpolate HRIRs ===============================================
% 1. Interpolate between close HRIRs
idx0 = 181; %index for 0° azimuth
idx1 = 182; %index for 1° azimuth
idx2 = 183; %index for 2° azimuth
x0_close = [hrtf.SourcePosition(idx0,:).' hrtf.SourcePosition(idx2,:).']/180*pi;
% weights (for target point: xs_close = hrtf.SourcePosition(idx1,:).'/180*pi )
weights_close = [.5 .5];
hrir_close = [hrtf.Data.IR(idx0,:,:); hrtf.Data.IR(idx2,:,:)];
hrir_close_ref = hrtf.Data.IR(idx1,:,:);

conf.ir.interpolationmethod = 'simple';
hrir_close_simple = interpolate_ir(hrir_close,weights_close,conf);
conf.ir.interpolationmethod = 'freqdomain';
hrir_close_fd = interpolate_ir(hrir_close,weights_close,conf);

% 2. Interpolate between distant HRIRs
idx0 = 181; %index for 0° azimuth
idx30 = 211; %index for 30° azimuth
idx60 = 241; %index for 60° azimuth
x0_dist = [hrtf.SourcePosition(idx0,:).' hrtf.SourcePosition(idx60,:).']/180*pi;
% weights (for target point: xs_dist = hrtf.SourcePosition(idx30,:).'/180*pi )
weights_dist = [.5 .5];
hrir_dist = [hrtf.Data.IR(idx0,:,:); hrtf.Data.IR(idx60,:,:)];
hrir_dist_ref = hrtf.Data.IR(idx30,:,:);

conf.ir.interpolationmethod = 'simple';
hrir_dist_simple = interpolate_ir(hrir_dist,weights_dist,conf);
conf.ir.interpolationmethod = 'freqdomain';
hrir_dist_fd = interpolate_ir(hrir_dist,weights_dist,conf);

% Plots
% impulse responses
figure
    plot(0:hrtf.API.N-1,squeeze(hrir_close(1,1,:)),'k'), hold on
    plot(0:hrtf.API.N-1,squeeze(hrir_close(2,1,:)),'b')
    plot(0:hrtf.API.N-1,squeeze(hrir_close_ref(1,1,:)),'g')
    plot(0:hrtf.API.N-1,squeeze(hrir_close_simple(1,1,:)),'r')
    plot(0:hrtf.API.N-1,squeeze(hrir_close_fd(1,1,:)),'m')
    grid
    xlabel('samples'), ylabel('amplitude')
    axis([0 160 -0.6 0.6])
    legend('hrir_1','hrir_2','hrir_{ref}','simple interp','freqdomain interp')
    title('Interpolation of close HRIRs')

figure
    plot(0:hrtf.API.N-1,squeeze(hrir_dist(1,1,:)),'k'), hold on
    plot(0:hrtf.API.N-1,squeeze(hrir_dist(2,1,:)),'b')
    plot(0:hrtf.API.N-1,squeeze(hrir_dist_ref(1,1,:)),'g')
    plot(0:hrtf.API.N-1,squeeze(hrir_dist_simple(1,1,:)),'r')
    plot(0:hrtf.API.N-1,squeeze(hrir_dist_fd(1,1,:)),'m')
    grid
    xlabel('samples'), ylabel('amplitude')
    axis([0 160 -0.6 0.6])
    legend('hrir_1','hrir_2','hrir_{ref}','simple interp','freqdomain interp')
    title('Interpolation of distant HRIRs')

% magnitude responses
f = (0:hrtf.API.N-1)/hrtf.API.N*hrtf.Data.SamplingRate;
figure
    semilogx(f,db(abs(fft(squeeze(hrir_close(1,1,:))))),'k'), hold on
    semilogx(f,db(abs(fft(squeeze(hrir_close(2,1,:))))),'b')
    semilogx(f,db(abs(fft(squeeze(hrir_close_ref(1,1,:))))),'g')
    semilogx(f,db(abs(fft(squeeze(hrir_close_simple(1,1,:))))),'r')
    semilogx(f,db(abs(fft(squeeze(hrir_close_fd(1,1,:))))),'m')
    grid
    xlabel('frequency in Hz'), ylabel('amplitude in dB')
    axis([0 hrtf.Data.SamplingRate/2 -60 20])
    legend('hrir_1','hrir_2','hrir_{ref}','simple interp','freqdomain interp',...
        'Location','SW')
    title('Interpolation of close HRIRs')

figure
    semilogx(f,db(abs(fft(squeeze(hrir_dist(1,1,:))))),'k'), hold on
    semilogx(f,db(abs(fft(squeeze(hrir_dist(2,1,:))))),'b')
    semilogx(f,db(abs(fft(squeeze(hrir_dist_ref(1,1,:))))),'g')
    semilogx(f,db(abs(fft(squeeze(hrir_dist_simple(1,1,:))))),'r')
    semilogx(f,db(abs(fft(squeeze(hrir_dist_fd(1,1,:))))),'m')
    grid
    xlabel('frequency in Hz'), ylabel('amplitude in dB')
    axis([0 hrtf.Data.SamplingRate/2 -60 20])
    legend('hrir_1','hrir_2','hrir_{ref}','simple interp','freqdomain interp',...
        'Location','SW')
    title('Interpolation of distant HRIRs')

status = true;
