%ANDROID_BRS_SET
%
%   Make a BRS set of shortend IRs for the Android phone
%

% AUTHOR: Hagen Wierstorf


%% ===== Variables ======================================================

% IR set to use
irset = 'Wittek_Studio_echoic';
%irset = 'FABIAN_postprocessed_anechoic';
% Target sampling rate
fs = 22050;
% Target length
nsamples = 128;
% Target angles
angles = (360:-1:1)/180*pi;

% File to save the BRS set
outfile = 'D:\data\measurements\BRIRs\android_128_Wittek_Studio_echoic_brs.wav';
%outfile = 'D:\data\measurements\HRIRs\android_128_FABIAN_pinta_anechoic_brs.wav';


%% ===== Configuration ==================================================

conf = SFS_config;

% Explicit configuration to ensure reproducibility
% === General audio ===
% Samplingrate (Hz)
conf.fs = 44100;
% Speed of sound (m/s)
conf.c = 343;


% === Plotting ===
% Plot the results etc.
conf.useplot = false;
% Use gnuplot
conf.usegnuplot = false;


%% ===== Computation ====================================================

irs = create_android_irs_mat(irset,nsamples,fs,conf);
brs = zeros(nsamples,2*length(angles));

% Generate a set of short IR
for ii = 1:length(angles)
    brs(:,(ii-1)*2+1:ii*2) = get_ir(irs,angles(ii)); 
end

% Normalize the IR set
brs = brs ./ (max(abs(brs(:)))+0.01);


%% ===== Save results ===================================================

% Write HRIRs to file
wavwrite(brs,conf.fs,16,outfile);
