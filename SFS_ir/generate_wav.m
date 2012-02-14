% Matlab script to generate wav files suitable for the SSR of the irs filesin this
% directory
%
%  For the SoundScape renderer (SSR) have a look at: http://tu-berlin.de/?id=ssr
%
%  NOTE: at the moment the script works only if your HRIR data were measured in
%  the whole hemisphere with equiangular resolution

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


% get files
files = dir('*.mat');

for ii = 1:length(files)
    irs = read_irs(files(ii).name);
    [len,ch] = size(irs.left);
    wav = zeros(len,2*ch);
    resolution = 360/ch;
    for jj = 1:ch
        ir = get_ir(irs,rad(jj*resolution));
        wav(:,2*jj-1) = ir(:,1);
        wav(:,2*jj)   = ir(:,2);
    end
    [directory, name, ext] = fileparts(files(ii).name);
    % Save as 32-bit float wave
    wavwrite(wav,irs.fs,32,[name,'.wav']);
end
