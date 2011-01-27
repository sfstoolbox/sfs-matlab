% Matlab script to generate an openDAFF file for the measured HRIR in this
% directory
%
%  For OpenDAFF see http://www.opendaff.org and the DAGA paper from Frank Wefers
%

% AUTHOR: Hagen Wierstorf

%% ===== Create custom metadata ==========================================
metadata = [];
metadata = daff_metadata_addKey(metadata, 'Label_Channel_1', 'String', ...
    'HRIR left ear');
metadata = daff_metadata_addKey(metadata, 'Label_Channel_2', 'String', ...
    'HRIR right ear');
metadata = daff_metadata_addKey(metadata, 'Manikin', 'String', ...
    'KEMAR Type 45BA');
metadata = daff_metadata_addKey(metadata, 'Loudspeaker', 'String', ...
    'Genelec 8030A');
metadata = daff_metadata_addKey(metadata, 'Room', 'String', ...
    'RAR, ITA, TU Berlin, Einsteinufer 13');
metadata = daff_metadata_addKey(metadata, 'Ears', 'String', ...
    'KEMAR large ears');
metadata = daff_metadata_addKey(metadata, 'Rotation', 'String', ...
    'Torso');
metadata = daff_metadata_addKey(metadata, 'Rotation_Hardware', 'String', ...
    'VariSphere');
metadata = daff_metadata_addKey(metadata, 'Exciter', 'String', ...
    'Emphasis_FFT18_44.1K');
metadata = daff_metadata_addKey(metadata, 'Excitation_Time', 'Float', ...
    5.3043);
metadata = daff_metadata_addKey(metadata, 'Date', 'String', ...
    '17-Nov-2010');


%% ===== Data function to extract the data from the mat files ============
function [data,samplerate,metadata] = daff_datafunc(alpha,beta,basepath)
    % Shift the angle
    alpha = -alpha+180;
    if alpha<=0
        alpha = alpha+360;
    end
    filestr = sprintf('%s/KEMAR_1deg_0.5m_large_ears%03.0f_%i.mat',...
        basepath,alpha+1,alpha);
    load(filestr);
    % We have only data in the horizontal plane
    if beta==90
        data = [vspolardata.ir_ch1(1:2048); vspolardata.ir_ch2(1:2048)];
    else
        % Dirac
        data = zeros(2,2048);
        data(:,1) = 1;
    end
    samplerate = 44100;
    metadata = [];
end


%% ===== Write the DAFF file =============================================
% Get dir of impulse responses
basedir = pwd;
daff_write('filename', 'KEMAR_1deg_0.5m_large_ears.daff', ...
    'content', 'ir', ...
    'datafunc', @daff_datafunc, ...
    'orient', [0 0 0], ...
    'basepath', basedir, ...
    'metadata', metadata, ...
    'mdist', 3.0, ... % distance between KEMAR and loudspeaker
    'samplerate', 44100, ...
    'quantization', 'float32', ... % NOTE: int16 or int32 for Octave only
    'channels', 2, ...
    'zthreshold', -200, ...  % Replace data with zeros for smaller amplitudes
    'alphares', 1, ...
    'betares', 90);


