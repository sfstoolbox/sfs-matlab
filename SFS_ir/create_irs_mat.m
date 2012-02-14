function irs = create_irs_mat(irset)
%CREATE_IRS_MAT Save a IR dataset as a mat file
%   Usage: create_irs_mat(irset)
%
%   Input parameters:
%       irset - Name of the IR set to be stored as a mat file. The name is used
%               in this file to identify the code to run and is used as a name
%               for the output mat file containing the IR set.
%
%   Output paramteres:
%       irs   - struct containing
%           .left           - left ear IR signal
%           .right          - right ear IR signal
%           .angle          - [azimuth; elevation] angles (rad)
%           .r0             - measurement distance of IR dataset (m)
%           .tag            - 'HRIR' or 'BRIR'
%           .description    - description of the IR data set
%
%   CREATE_IRS_MAT(irset) reads a IR dataset from wav file etc. and stores it in
%   a mat file in the format needed by SFS. If you have your own IR dataset you
%   will add a section to this file in order to read it and store it in the
%   desired format.
%
%   IR sets, that are supported at the moment are (given by irset name):
%
%       'FABIAN_anechoic'
%           HRIR data set measured by Alexander Lindau in the
%           "Reflexionsarmen Raum im ITA TU Berlin" using the FABIAN
%           manikin. The set contains only data for the horizontal plane
%           with an azimuth resolution of 1.25°
%           NOTE: the real azimuth resolution of the FABIAN is 2.5°. 1.25° has
%           been reached by manual move FABIAN 1.25°. Therefore we have an
%           systematical error in the data.
%
%       'Wittek_Studio_echoic'
%           BRIR data set measured by Robert Wittek in Studio (Bochum) with
%           the ... manikin. The set contains only data for the horizontal
%           plane with an azimuth resolution of 1° (NOTE: if you need a
%           finer resolution have a look at the original data, there Wittek
%           uses a resolution of just 0.1°)
%
%
%   see also: read_irs, get_ir, ir_intpol
%

% AUTHOR: Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Configuration ==================================================
% Path to the QU files
QU_PATH = '/home/hagen/data';


%% ===== Checking of input  parameters ==================================

if nargchk(1,1,nargin)
    error('Wrong number of args. Usage: irs = create_irs_mat(irset)');
end
if ~ischar(irset)
    error('%s: irset has to be a string!',upper(mfilename));
end


%% ===== Read HRTF files ================================================

% === Alex Lindau - Pinta HRIRs ===
if strcmp(irset,'FABIAN_postprocessed_anechoic')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Alex Lindaus anechoic measurements with FABIAN in Pinta. ',...
         'Used elevation angle: 0°; azimuth resolution: 1.25°.'];
    irs.tag = 'HRIR';
    % Measurement distance
    irs.r0 = 2.5;
    % Configuration for the given HRIR set
    srcdir = sprintf(...
        '%s/measurements/HRIRs/FABIAN_anechoic/postprocessed/source1/',...
        QU_PATH);
    outdir = sprintf('%s/measurements/HRIRs',QU_PATH);
    % Angular stepsize
    step = 1.25;

    % Read the data
    count=0;
    for phi = -180:step:180-step
        count=count+1;

        if phi<0
            src=strcat(srcdir,'N');
        else
            src=strcat(srcdir,'P');
        end

        if phi<=-100 || phi>=100
            num=strcat(num2str(abs(phi*100)));
        end

        if phi>-100 && phi<100
            pad='0';
            num=strcat(pad,num2str(abs(phi*100)));
        end

        if phi>-10 && phi<10
            pad='00';
            num=strcat(pad,num2str(abs(phi*100)));
        end

        if phi>-1 && phi<1
            pad='0000';
            num=strcat(pad,num2str(abs(phi*100)));
        end

        ir=wavread(strcat(src,num,'.wav'));

        irs.angle(:,count) = [phi / 180*pi; 0];
        irs.left(:,count) = ir(:,1);
        irs.right(:,count) = ir(:,2);

    end

    % Reorder the HRTF set
    % For HRTF interpolation we need an increasing angle order in the HRTF set,
    % so we have to reorder them
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);

elseif strcmp(irset,'FABIAN_anechoic')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Alex Lindaus anechoic measurements with FABIAN in Pinta. ',...
         'Used elevation angle: 0°; azimuth resolution: 1.25°.'];
    irs.tag = 'HRIR';
    % Measurement distance
    irs.r0 = 2.5;
    % Configuration for the given HRIR set
    srcdir = sprintf(...
        '%s/measurements/HRIRs/FABIAN_anechoic/original/source1/',QU_PATH);
    outdir = sprintf('%s/measurements/HRIRs',QU_PATH);
    % Angular stepsize
    step = 1.25;

    % Read the data
    count=0;
    for phi = -180:step:180-step
        count=count+1;

        if phi<0
            src=strcat(srcdir,'N');
        else
            src=strcat(srcdir,'P');
        end

        if phi<=-100 || phi>=100
            num=strcat(num2str(abs(phi*100)));
        end

        if phi>-100 && phi<100
            pad='0';
            num=strcat(pad,num2str(abs(phi*100)));
        end

        if phi>-10 && phi<10
            pad='00';
            num=strcat(pad,num2str(abs(phi*100)));
        end

        if phi>-1 && phi<1
            pad='0000';
            num=strcat(pad,num2str(abs(phi*100)));
        end

        ir=wavread(strcat(src,num,'.wav'));

        irs.angle(:,count) = [phi / 180*pi; 0];
        irs.left(:,count) = ir(:,1);
        irs.right(:,count) = ir(:,2);

    end

    % Reorder the HRTF set
    % For HRTF interpolation we need an increasing angle order in the HRTF set,
    % so we have to reorder them
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);


% === Robert Wittek - BRIRs ===
elseif strcmp(irset,'Wittek_Studio_echoic')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Robert Witteks measurements with ... in Studio. ',...
         'Used elevation angle: 0°; azimuth resolution: 1°.'];
    irs.tag = 'BRIR';
    % Measurement distance
    irs.r0 = 1;
    % Configuration for the given HRIR set
    srcdir = sprintf('%s/measurements/BRIRs/Wittek',QU_PATH);
    outdir = sprintf('%s/measurements/BRIRs',QU_PATH);
    % Angular stepsize
    step = 1;

    % Disable warnings, because the hrtf struct of Wittek has the class
    % :all: that is no longer supported (?)
    warning('off', 'all')

    % Read the data
    count=1;

    % NOTE: in order to have an ascending order of the angles in the IR
    % data set, we are using two different loops to read the data, starting
    % at -pi and ending at pi-eps

    % -pi..0-eps
    for phi = 1800:step*10:3590

        irfile = sprintf('%s/Studio360_KHneu_L3600_%i.mat',srcdir,phi);
        load(irfile);
        irs.angle(:,count) = [correct_angle(phi/1800*pi); 0];
        irs.left(:,count) = hrtf.data(:,1);
        irs.right(:,count) = hrtf.data(:,2);
        count = count+1;

    end

    % 0..pi-eps
    for phi = 0:step*10:1790

        irfile = sprintf('%s/Studio360_KHneu_L3600_%i.mat',srcdir,phi);
        load(irfile);
        irs.angle(:,count) = [phi/1800*pi; 0];
        irs.left(:,count) = hrtf.data(:,1);
        irs.right(:,count) = hrtf.data(:,2);
        count = count+1;

    end

    % Reenable the warnings
    warning('on', 'all')

% === KEMAR RAR 1deg 3m small ===
elseif strcmp(irset,'KEMAR_RAR_1deg_3m_small')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Measurements with KEMAR using small ears in RAR of TU Berlin. ',...
         'Used elevation angle: 0°; azimuth resolution: 1°.'];
    irs.tag = 'HRIR';
    % Measurement distance (m)
    irs.r0 = 3;
    % Configuration for the given HRIR set
    srcdir = sprintf(...
        '%s/measurements/HRIRs/VS_POLAR/KEMAR_1deg_3m_small_ears',QU_PATH);
    outdir = sprintf('%s/measurements/HRIRs',QU_PATH);
    % Angular stepsize
    step = 1;

    for ii=1:360
        irfile = sprintf('%s/KEMAR_1deg_3m_%03.0f_%i.mat',srcdir,ii,ii-1);
        load(irfile);
        irs.angle(:,ii) = [correct_angle((-ii-1+180)/180*pi); 0];
        irs.left(:,ii) = vspolardata.ir_ch1;
        irs.right(:,ii) = vspolardata.ir_ch2;
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);


% === KEMAR RAR 1deg 3m control ===
elseif strcmp(irset,'KEMAR_RAR_1deg_3m_small_control')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Control measurements with KEMAR using small ears', ...
         'in RAR of TU Berlin. ',...
         'Used elevation angle: 0°; azimuth resolution: 1°.'];
    irs.tag = 'HRIR';
    % Measurement distance (m)
    irs.r0 = 3;
    % Configuration for the given HRIR set
    srcdir = sprintf(...
        '%s/measurements/HRIRs/VS_POLAR/KEMAR_1deg_3m_small_ears_control_',...
        QU_PATH);
    outdir = sprintf('%s/measurements/HRIRs',QU_PATH);
    % Angular stepsize
    step = 1;

    for ii=1:360
        irfile = sprintf('%s/KEMAR_1deg_3m_control_%03.0f_%i.mat',srcdir,ii,ii-1);
        load(irfile);
        irs.angle(:,ii) = [correct_angle((-ii-1+180)/180*pi); 0];
        irs.left(:,ii) = vspolardata.ir_ch1;
        irs.right(:,ii) = vspolardata.ir_ch2;
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);

% === KEMAR RAR 1deg 3m large ear ===
elseif strcmp(irset,'KEMAR_RAR_1deg_3m_large')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Measurements with KEMAR using large ears in RAR of TU Berlin. ',...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg.'];
    irs.tag = 'HRIR';
    % Measurement distance (m)
    irs.r0 = 3;
    % Configuration for the given HRIR set
    srcdir = sprintf(...
        '%s/measurements/HRIRs/KEMAR_1deg_3m_large_ears_',QU_PATH);
    outdir = sprintf('%s/measurements/HRIRs',QU_PATH);
    % Angular stepsize
    step = 1;

    for ii=1:360
        irfile = sprintf('%s/KEMAR_1deg_3m_large_ears_%03.0f_%i.mat',...
            srcdir,ii,ii-1);
        load(irfile);
        irs.angle(:,ii) = [correct_angle((-ii-1+180)/180*pi); 0];
        irs.left(:,ii) = vspolardata.ir_ch1;
        irs.right(:,ii) = vspolardata.ir_ch2;
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);

% === KEMAR RAR 1deg 2m large ear ===
elseif strcmp(irset,'KEMAR_RAR_1deg_2m_large')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Measurements with KEMAR using large ears in RAR of TU Berlin. ',...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg.'];
    irs.tag = 'HRIR';
    % Measurement distance (m)
    irs.r0 = 3;
    % Configuration for the given HRIR set
    srcdir = sprintf(...
        '%s/measurements/HRIRs/KEMAR_1deg_2m_large_ears_',QU_PATH);
    outdir = sprintf('%s/measurements/HRIRs',QU_PATH);
    % Angular stepsize
    step = 1;

    for ii=1:360
        irfile = sprintf('%s/KEMAR_1deg_2m_large_ears_%03.0f_%i.mat',...
            srcdir,ii,ii-1);
        load(irfile);
        irs.angle(:,ii) = [correct_angle((-ii-1+180)/180*pi); 0];
        irs.left(:,ii) = vspolardata.ir_ch1;
        irs.right(:,ii) = vspolardata.ir_ch2;
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);

% === KEMAR RAR 1deg 3m large ear ===
elseif strcmp(irset,'KEMAR_RAR_1deg_1m_large')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Measurements with KEMAR using large ears in RAR of TU Berlin. ',...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg.'];
    irs.tag = 'HRIR';
    % Measurement distance (m)
    irs.r0 = 3;
    % Configuration for the given HRIR set
    srcdir = sprintf(...
        '%s/measurements/HRIRs/VS_POLAR/KEMAR_1deg_1m_large_ear',QU_PATH);
    outdir = sprintf('%s/measurements/HRIRs',QU_PATH);
    % Angular stepsize
    step = 1;

    for ii=1:360
        irfile = sprintf('%s/KEMAR_1deg_1m_large_ear%03.0f_%i.mat',...
            srcdir,ii,ii-1);
        load(irfile);
        irs.angle(:,ii) = [correct_angle((-ii-1+180)/180*pi); 0];
        irs.left(:,ii) = vspolardata.ir_ch1;
        irs.right(:,ii) = vspolardata.ir_ch2;
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);

% === KEMAR RAR 1deg 0.5m large ear ===
elseif strcmp(irset,'KEMAR_RAR_1deg_0.5m_large')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Measurements with KEMAR using large ears in RAR of TU Berlin. ',...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg.'];
    irs.tag = 'HRIR';
    % Measurement distance (m)
    irs.r0 = 3;
    % Configuration for the given HRIR set
    srcdir = sprintf(...
        '%s/measurements/HRIRs/VS_POLAR/KEMAR_1deg_0.5m_large_ears',QU_PATH);
    outdir = sprintf('%s/measurements/HRIRs',QU_PATH);
    % Angular stepsize
    step = 1;

    for ii=1:360
        irfile = sprintf('%s/KEMAR_1deg_0.5m_large_ears%03.0f_%i.mat',...
            srcdir,ii,ii-1);
        load(irfile);
        irs.angle(:,ii) = [correct_angle((-ii-1+180)/180*pi); 0];
        irs.left(:,ii) = vspolardata.ir_ch1;
        irs.right(:,ii) = vspolardata.ir_ch2;
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);

% === AUDIS HRIR ===
elseif strcmp(irset,'AUDIS')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['AUDIS HRIR set. ',...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg.'];
    irs.tag = 'HRIR';
    % Measurement distance (m)
    irs.r0 = NaN;
    % Configuration for the given HRIR set
    srcdir = sprintf(...
        '%s/measurements/HRIRs/AUDIS/HRIRs_AUDIS',QU_PATH);
    outdir = sprintf('%s/measurements/HRIRs',QU_PATH);
    % Angular stepsize
    step = 1;

    for ii=1:360
        irfile = sprintf('%s/hrir%03.0f.wav',srcdir,ii-1);
        ir = wavread(irfile);
        irs.angle(:,ii) = [correct_angle(ii/180*pi); 0];
        irs.left(:,ii) = ir(:,1);
        irs.right(:,ii) = ir(:,2);
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);

% === Aachen BRIR ===
elseif strcmp(irset,'Aachen')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Aachen BRIR set for Stairway. ',...
         'Used elevation angle: 0 deg; azimuth resolution: 15 deg.'];
    irs.tag = 'BRIR';
    outdir = sprintf('%s/measurements/BRIRs',QU_PATH);
    addpath([QU_PATH ...
        '/measurements/BRIRs/Aachen_Impulse_Response_Database_v1.2']);

    % === Stairway ===
    % Measurement distance (m)
    irs.r0 = 3;

    airpar.fs = 44100;      % samplingrate
    airpar.rir_type = 1;    % BRIR
    airpar.room = 5;        % Stairway
    airpar.head = 1;        % with dummy head
    airpar.rir_no = 3;      % r = 3m


    % Angular stepsize
    step = 15;

    count = 1;
    for ii=15:15:180
        airpar.azimuth = ii;
        irs.angle(:,count) = [correct_angle((-ii+90)/180*pi); 0];
        airpar.channel = 1;
        irs.left(:,count) = load_air(airpar)';
        airpar.channel = 0;
        irs.right(:,count) = load_air(airpar)';
        count = count+1;
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);

% === MIT KEMAR ===
elseif strcmp(irset,'KEMAR_MIT')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Measurements with KEMAR at MIT. ',...
         'See: http://sound.media.mit.edu/resources/KEMAR.html',...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg.'];
    irs.tag = 'HRIR';
    % Measurement distance (m)
    irs.r0 = 1.7;
    % Configuration for the given HRIR set
    srcdir = sprintf(...
        '%s/measurements/HRIRs/MIT_Kemar/compact/elev0',QU_PATH);
    outdir = sprintf('%s/measurements/HRIRs',QU_PATH);
    % Angular stepsize
    step = 5;

    count = 1;
    for ii=0:5:180
        irfile = sprintf('%s/H0e%03.0fa.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle(-ii/180*pi); 0];
        irs.left(:,count) = ir(:,1);
        irs.right(:,count) = ir(:,2);
        count = count+1;
    end
    for ii=5:5:175
        irfile = sprintf('%s/H0e%03.0fa.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle(ii/180*pi); 0];
        irs.left(:,count) = ir(:,2);
        irs.right(:,count) = ir(:,1);
        count = count+1;
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);


% === FABIAN Sputnik ===
elseif strcmp(irset,'SPUTNIK1')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Measurements with FABIAN at Sputnik. ',...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg.'];
    irs.tag = 'BRIR';
    % Measurement distance (m)
    irs.r0 = NaN;
    % Configuration for the given BRIR set
    srcdir = sprintf(...
        '%s/measurements/BRIRs/Sputnik_corrected/source1',QU_PATH);
    outdir = sprintf('%s/measurements/BRIRs',QU_PATH);
    % Angular stepsize
    step = 1;

    count = 1;
    ir_length = 44100;
    for ii=1:1:90
        irfile = sprintf('%s/P00_N%02.0f0.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle(-(ii-91)/180*pi); 0];
        irs.left(:,count) = ir(1:ir_length,1);
        irs.right(:,count) = ir(1:ir_length,2);
        count = count+1;
    end
    for ii=0:1:90
        irfile = sprintf('%s/P00_P%02.0f0.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle((ii-90)/180*pi); 0];
        irs.left(:,count) = ir(1:ir_length,2);
        irs.right(:,count) = ir(1:ir_length,1);
        count = count+1;
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);

% === FABIAN Sputnik ===
elseif strcmp(irset,'SPUTNIK2')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Measurements with FABIAN at Sputnik. ',...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg.'];
    irs.tag = 'BRIR';
    % Measurement distance (m)
    irs.r0 = NaN;
    % Configuration for the given BRIR set
    srcdir = sprintf(...
        '%s/measurements/BRIRs/Sputnik_corrected/source2',QU_PATH);
    outdir = sprintf('%s/measurements/BRIRs',QU_PATH);
    % Angular stepsize
    step = 1;

    count = 1;
    ir_length = 44100;
    for ii=1:1:90
        irfile = sprintf('%s/P00_N%02.0f0.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle(-(ii-91)/180*pi); 0];
        irs.left(:,count) = ir(1:ir_length,1);
        irs.right(:,count) = ir(1:ir_length,2);
        count = count+1;
    end
    for ii=0:1:90
        irfile = sprintf('%s/P00_P%02.0f0.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle((ii-90)/180*pi); 0];
        irs.left(:,count) = ir(1:ir_length,2);
        irs.right(:,count) = ir(1:ir_length,1);
        count = count+1;
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);

% === FABIAN Sputnik ===
elseif strcmp(irset,'SPUTNIK3')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Measurements with FABIAN at Sputnik. ',...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg.'];
    irs.tag = 'BRIR';
    % Measurement distance (m)
    irs.r0 = NaN;
    % Configuration for the given BRIR set
    srcdir = sprintf(...
        '%s/measurements/BRIRs/Sputnik_corrected/source3',QU_PATH);
    outdir = sprintf('%s/measurements/BRIRs',QU_PATH);
    % Angular stepsize
    step = 1;

    count = 1;
    ir_length = 44100;
    for ii=1:1:90
        irfile = sprintf('%s/P00_N%02.0f0.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle(-(ii-91)/180*pi); 0];
        irs.left(:,count) = ir(1:ir_length,1);
        irs.right(:,count) = ir(1:ir_length,2);
        count = count+1;
    end
    for ii=0:1:90
        irfile = sprintf('%s/P00_P%02.0f0.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle((ii-90)/180*pi); 0];
        irs.left(:,count) = ir(1:ir_length,2);
        irs.right(:,count) = ir(1:ir_length,1);
        count = count+1;
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);

% === FABIAN Sputnik ===
elseif strcmp(irset,'SPUTNIK4')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Measurements with FABIAN at Sputnik. ',...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg.'];
    irs.tag = 'BRIR';
    % Measurement distance (m)
    irs.r0 = NaN;
    % Configuration for the given BRIR set
    srcdir = sprintf(...
        '%s/measurements/BRIRs/Sputnik_corrected/source4',QU_PATH);
    outdir = sprintf('%s/measurements/BRIRs',QU_PATH);
    % Angular stepsize
    step = 1;

    count = 1;
    ir_length = 44100;
    for ii=1:1:90
        irfile = sprintf('%s/P00_N%02.0f0.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle(-(ii-91)/180*pi); 0];
        irs.left(:,count) = ir(1:ir_length,1);
        irs.right(:,count) = ir(1:ir_length,2);
        count = count+1;
    end
    for ii=0:1:90
        irfile = sprintf('%s/P00_P%02.0f0.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle((ii-90)/180*pi); 0];
        irs.left(:,count) = ir(1:ir_length,2);
        irs.right(:,count) = ir(1:ir_length,1);
        count = count+1;
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);

% === FABIAN Sputnik ===
elseif strcmp(irset,'SPUTNIK5')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Measurements with FABIAN at Sputnik. ',...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg.'];
    irs.tag = 'BRIR';
    % Measurement distance (m)
    irs.r0 = NaN;
    % Configuration for the given BRIR set
    srcdir = sprintf(...
        '%s/measurements/BRIRs/Sputnik_corrected/source5',QU_PATH);
    outdir = sprintf('%s/measurements/BRIRs',QU_PATH);
    % Angular stepsize
    step = 1;

    count = 1;
    ir_length = 44100;
    for ii=1:1:90
        irfile = sprintf('%s/P00_N%02.0f0.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle(-(ii-91)/180*pi); 0];
        irs.left(:,count) = ir(1:ir_length,1);
        irs.right(:,count) = ir(1:ir_length,2);
        count = count+1;
    end
    for ii=0:1:90
        irfile = sprintf('%s/P00_P%02.0f0.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle((ii-90)/180*pi); 0];
        irs.left(:,count) = ir(1:ir_length,2);
        irs.right(:,count) = ir(1:ir_length,1);
        count = count+1;
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);

% === FABIAN Sputnik ===
elseif strcmp(irset,'SPUTNIK6')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Measurements with FABIAN at Sputnik. ',...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg.'];
    irs.tag = 'BRIR';
    % Measurement distance (m)
    irs.r0 = NaN;
    % Configuration for the given BRIR set
    srcdir = sprintf(...
        '%s/measurements/BRIRs/Sputnik_corrected/source6',QU_PATH);
    outdir = sprintf('%s/measurements/BRIRs',QU_PATH);
    % Angular stepsize
    step = 1;

    count = 1;
    ir_length = 44100;
    for ii=1:1:90
        irfile = sprintf('%s/P00_N%02.0f0.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle(-(ii-91)/180*pi); 0];
        irs.left(:,count) = ir(1:ir_length,1);
        irs.right(:,count) = ir(1:ir_length,2);
        count = count+1;
    end
    for ii=0:1:90
        irfile = sprintf('%s/P00_P%02.0f0.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle((ii-90)/180*pi); 0];
        irs.left(:,count) = ir(1:ir_length,2);
        irs.right(:,count) = ir(1:ir_length,1);
        count = count+1;
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);

% === FABIAN Sputnik ===
elseif strcmp(irset,'SPUTNIK7')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Measurements with FABIAN at Sputnik. ',...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg.'];
    irs.tag = 'BRIR';
    % Measurement distance (m)
    irs.r0 = NaN;
    % Configuration for the given BRIR set
    srcdir = sprintf(...
        '%s/measurements/BRIRs/Sputnik_corrected/source7',QU_PATH);
    outdir = sprintf('%s/measurements/BRIRs',QU_PATH);
    % Angular stepsize
    step = 1;

    count = 1;
    ir_length = 44100;
    for ii=1:1:90
        irfile = sprintf('%s/P00_N%02.0f0.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle(-(ii-91)/180*pi); 0];
        irs.left(:,count) = ir(1:ir_length,1);
        irs.right(:,count) = ir(1:ir_length,2);
        count = count+1;
    end
    for ii=0:1:90
        irfile = sprintf('%s/P00_P%02.0f0.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle((ii-90)/180*pi); 0];
        irs.left(:,count) = ir(1:ir_length,2);
        irs.right(:,count) = ir(1:ir_length,1);
        count = count+1;
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);

% === FABIAN Sputnik ===
elseif strcmp(irset,'SPUTNIK8')

    % Initialize a new IR struct
    irs = {};
    % Description
    irs.description = ...
        ['Measurements with FABIAN at Sputnik. ',...
         'Used elevation angle: 0 deg; azimuth resolution: 1 deg.'];
    irs.tag = 'BRIR';
    % Measurement distance (m)
    irs.r0 = NaN;
    % Configuration for the given BRIR set
    srcdir = sprintf(...
        '%s/measurements/BRIRs/Sputnik_corrected/source8',QU_PATH);
    outdir = sprintf('%s/measurements/BRIRs',QU_PATH);
    % Angular stepsize
    step = 1;

    count = 1;
    ir_length = 44100;
    for ii=1:1:90
        irfile = sprintf('%s/P00_N%02.0f0.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle(-(ii-91)/180*pi); 0];
        irs.left(:,count) = ir(1:ir_length,1);
        irs.right(:,count) = ir(1:ir_length,2);
        count = count+1;
    end
    for ii=0:1:89
        irfile = sprintf('%s/P00_P%02.0f0.wav',srcdir,ii);
        ir = wavread(irfile);
        irs.angle(:,count) = [correct_angle((ii-89)/180*pi); 0];
        irs.left(:,count) = ir(1:ir_length,2);
        irs.right(:,count) = ir(1:ir_length,1);
        count = count+1;
    end

    % Resort the data
    hphi = irs.angle(1,:);
    [hphi,idx] = sort(hphi);
    irs.angle = irs.angle(:,idx);
    irs.left = irs.left(:,idx);
    irs.right = irs.right(:,idx);

else
    error('The given IR set is unknown');
end


%% ===== Write the HRIR mat file =========================================
outfile = sprintf('%s/%s.mat',outdir,irset);
save('-v7',outfile,'irs');

