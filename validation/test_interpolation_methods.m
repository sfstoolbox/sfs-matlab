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

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2019 SFS Toolbox Developers                             *
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
% https://sfs.readthedocs.io                            sfstoolbox@gmail.com *
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
weights = [.5; .5];
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
conf.ir.interpolationmethod = 'timedomain';
h_int1_td = interpolate_ir(irs1,weights,conf);

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
conf.ir.interpolationmethod = 'timedomain';
h_int2_td = interpolate_ir(irs2,weights,conf);

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
conf.ir.interpolationmethod = 'timedomain';
h_int3_td = interpolate_ir(irs3,weights,conf);

% Plots
% impulse responses
ir_list = {
    irs1(1,1,:)         , irs2(1,1,:)          , irs3(1,1,:)
    irs1(2,1,:)         , irs2(2,1,:)          , irs3(2,1,:)
    h_int1_simple(1,1,:), h_int2_simple(1,1,:) , h_int3_simple(1,1,:)
    h_int1_fd(1,1,:)    , h_int2_fd(1,1,:)     , h_int3_fd(1,1,:)
    h_int1_td(1,1,:)    , h_int2_td(1,1,:)     , h_int3_td(1,1,:)
    };  % list of IRs

col = {'k' 'b' 'r' 'm' 'c'};

titlelabels = {
    'Interpolation between Dirac impulses with one sample in between'
    'Interpolation of neighbouring Dirac impulses at impulse response start'
    'Interpolation of neighbouring Dirac impulses in middle of impulse response'
    };

t = 0:N-1;  % time-axis in samples

for cdx=1:size(ir_list,2)
    figure
    for rdx=1:size(ir_list,1)
        plot(t,squeeze(ir_list{rdx, cdx}), col{rdx})
        hold on
    end
    hold off
    grid
    xlabel('samples')
    ylabel('amplitude')
    legend('h_1','h_2','simple interp','freqdomain interp','timedomain interp')
    title(titlelabels{cdx})
end

%% ===== Interpolate HRIRs ===============================================
% 1. Interpolate between close HRIRs
idx0 = 181; %index for 0° azimuth
idx1 = 182; %index for 1° azimuth
idx2 = 183; %index for 2° azimuth
x0_close = [hrtf.SourcePosition(idx0,:).' hrtf.SourcePosition(idx2,:).']/180*pi;
% weights (for target point: xs_close = hrtf.SourcePosition(idx1,:).'/180*pi )
weights_close = [.5; .5];
hrir_close = [hrtf.Data.IR(idx0,:,:); hrtf.Data.IR(idx2,:,:)];
hrir_close_ref = hrtf.Data.IR(idx1,:,:);

conf.ir.interpolationmethod = 'simple';
hrir_close_simple = interpolate_ir(hrir_close,weights_close,conf);
conf.ir.interpolationmethod = 'freqdomain';
hrir_close_fd = interpolate_ir(hrir_close,weights_close,conf);
conf.ir.interpolationmethod = 'timedomain';
hrir_close_td = interpolate_ir(hrir_close,weights_close,conf);

% 2. Interpolate between distant HRIRs
idx0 = 181; %index for 0° azimuth
idx30 = 211; %index for 30° azimuth
idx60 = 241; %index for 60° azimuth
x0_dist = [hrtf.SourcePosition(idx0,:).' hrtf.SourcePosition(idx60,:).']/180*pi;
% weights (for target point: xs_dist = hrtf.SourcePosition(idx30,:).'/180*pi )
weights_dist = [.5; .5];
hrir_dist = [hrtf.Data.IR(idx0,:,:); hrtf.Data.IR(idx60,:,:)];
hrir_dist_ref = hrtf.Data.IR(idx30,:,:);

conf.ir.interpolationmethod = 'simple';
hrir_dist_simple = interpolate_ir(hrir_dist,weights_dist,conf);
conf.ir.interpolationmethod = 'freqdomain';
hrir_dist_fd = interpolate_ir(hrir_dist,weights_dist,conf);
conf.ir.interpolationmethod = 'timedomain';
hrir_dist_td = interpolate_ir(hrir_dist,weights_dist,conf);

% Plots
hrir_list = {
    hrir_close(1,1,:)       , hrir_dist(1,1,:)
    hrir_close(2,1,:)       , hrir_dist(2,1,:)
    hrir_close_ref(1,1,:)   , hrir_dist_ref(1,1,:)
    hrir_close_simple(1,1,:), hrir_dist_simple(1,1,:)
    hrir_close_fd(1,1,:)    , hrir_dist_fd(1,1,:)
    hrir_close_td(1,1,:)    , hrir_dist_td(1,1,:)
    };  % list of HRIRs

col = {'k' 'b' 'g' 'r' 'm' 'c'};
titlelabels = {'close', 'distant'};
legendpos = {'northeast', 'southwest'};

t = 0:hrtf.API.N-1;  % time-axis in samples
f = (0:hrtf.API.N-1)/hrtf.API.N*hrtf.Data.SamplingRate;  % frequency axis

for pdx=1:2  % impulse responses or magnitude spectrum
    for cdx=1:size(hrir_list,2)
        figure
        for rdx=1:size(hrir_list,1)
            if pdx == 1
                % impulse responses
                plot(t,squeeze(hrir_list{rdx, cdx}),col{rdx})
            else
                % magnitude spectra
                mag_spectrum = db(abs(fft(squeeze(hrir_list{rdx, cdx}))));
                semilogx(f(2:end), mag_spectrum(2:end), col{rdx});
            end
            hold on
        end
        hold off

        if pdx == 1
            % impulse responses
            xlabel('samples')
            ylabel('amplitude')
            axis([0 160 -0.6 0.6])
        else
            % magnitude spectra
            xlabel('frequency in Hz')
            ylabel('amplitude in dB')
            axis([f(2) hrtf.Data.SamplingRate/2 -60 20])
        end
        legend('hrir_1','hrir_2','hrir_{ref}','simple interp', ...
            'freqdomain interp', 'timedomain interp','Location', legendpos{pdx})
        title(['Interpolation of ' titlelabels{cdx} ' HRIRs'])
    end
end

status = true;
