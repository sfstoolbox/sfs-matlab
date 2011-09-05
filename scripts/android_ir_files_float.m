%CREATE_ANDROID_IRS
%
%   CREATE_ANDROID_IRS creates a short version of the given irsset 
%   for the Android phone. Therefore the single IRs are stored as binary 
%   files using float32 in the corresponding directory.
%

% AUTHOR: Hagen Wierstorf, Sascha Spors


%% ===== Variables ======================================================

%irset = '~/data/ir_databases/Wittek_KEMAR/Wittek_KEMAR_studio_src1_0deg.mat';
irset = '~/data/ir_databases/QU_KEMAR/QU_KEMAR_anechoic_AKGK601_3m.mat';
fs = 22050;
nsamples = 256;
outdir = sprintf('~/data/measurements/BRIRs/android_%i_KEMAR_headphone_comp',nsamples);

irs = read_irs(irset);
angles = irs.apparent_azimuth;


%% ===== Configuration ==================================================
% Load default configuration values
conf = SFS_config;


%% ===== Computation ====================================================

irs = create_android_irs_mat(irs,nsamples,fs,conf);


%% ===== Save results ===================================================

% Write IRs to file
for ii = 1:length(angles)
    
    ir = get_ir(irs,angles(ii));
    
    % Left HRIR signal
    pstr = sprintf('%s/HRIR_left_%d.dat',outdir,round(angles(ii)/pi*180));
    fid=fopen(pstr,'w');
    %fprintf(fid,'%d\r\n',short_hrirs(1:nsamples,ii,1));
    fwrite(fid,ir(:,1),'double');
    fclose(fid);
    
    % Right HRIR signal
    pstr = sprintf('%s/HRIR_right_%d.dat',outdir,round(angles(ii)/pi*180));
    fid=fopen(pstr,'w');
    %fprintf(fid,'%d\r\n',short_hrirs(1:nsamples,ii,2));
    fwrite(fid,ir(:,2),'double');
    fclose(fid);
    
end

