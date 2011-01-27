%CREATE_ANDROID_IRS
%
%   CREATE_ANDROID_IRS creates a short version of the given irsset 
%   for the Android phone. Therefore the single IRs are stored as binary 
%   files using float32 in the corresponding directory.
%

% AUTHOR: Hagen Wierstorf, Sascha Spors


%% ===== Variables ======================================================

% IR set to use
irset = 'Wittek_Studio_echoic';
%irset = 'FABIAN_pinta_anechoic';
% Target sampling rate
fs = 22050;
% Target length
nsamples = 128;
% Target angles
angles = (-180:1:180)/180*pi;

% Outdir for the dat files
outdir = 'D:\data\measurements\BRIRs\android_128_Wittek_Studio_echoic';
%outdir = 'D:\data\measurements\HRIRs\android_196_FABIAN_pinta_anechoic';


%% ===== Configuration ==================================================

% Load default configuration values
conf = SFS_config;


%% ===== Computation ====================================================

irs = create_android_irs_mat(irset,nsamples,fs,conf);


%% ===== Save results ===================================================

% Write IRs to file
for ii = 1:length(angles)
    
    ir = get_ir(irs,angles(ii));
    
    % Left HRIR signal
    pstr = sprintf('%s/HRTF_left_%d.dat',outdir,round(angles(ii)/pi*180));
    fid=fopen(pstr,'w');
    %fprintf(fid,'%d\r\n',short_hrirs(1:nsamples,ii,1));
    fwrite(fid,ir(:,1),'float32');
    fclose(fid);
    
    % Right HRIR signal
    pstr = sprintf('%s/HRTF_right_%d.dat',outdir,round(angles(ii)/pi*180));
    fid=fopen(pstr,'w');
    %fprintf(fid,'%d\r\n',short_hrirs(1:nsamples,ii,2));
    fwrite(fid,ir(:,2),'float32');
    fclose(fid);
    
end

