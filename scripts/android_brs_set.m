%ANDROID_BRS_SET
%
%   Make a BRS set of shortend IRs for the Android phone
%

% AUTHOR: Hagen Wierstorf


%% ===== Variables ======================================================

% IR set to use
%irset = '~/data/ir_databases/Wittek_KEMAR/Wittek_KEMAR_studio_src1_0deg.mat';
irset = '~/data/ir_databases/QU_KEMAR/QU_KEMAR_anechoic_AKGK601_3m.mat';
% Target sampling rate
fs = 22050;
% Target length
nsamples = 128;
% File to save the BRS set
outfile = sprintf('~/data/measurements/BRIRs/android_%i_KEMAR_headphone_comp_brs.wav', ...
    nsamples);

irs = read_irs(irset);

%% ===== Configuration ==================================================

conf = SFS_config;

% Explicit configuration to ensure reproducibility
% === General audio ===
% Samplingrate (Hz)
conf.fs = 44100;
% Speed of sound (m/s)
conf.c = 343;
angles = rad(conf.brsangles);

% === Plotting ===
% Plot the results etc.
conf.useplot = false;
% Use gnuplot
conf.usegnuplot = false;


%% ===== Computation ====================================================

irs = create_android_irs_mat(irs,nsamples,fs,conf);
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
