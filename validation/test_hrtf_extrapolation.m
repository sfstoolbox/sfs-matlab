function boolean = test_hrtf_extrapolation(hrtf_set)
%TEST_HRTF_EXTRAPOLATION tests the HRTF extrapolation functions
%
%   Usage: boolean = test_hrtf_extrapolation(hrtf_set)
%
%   Input parameters:
%       hrtf_set - HRTFs for testing extrapolation. Can be
%                    'QU_KEMAR'
%                    'FABIAN_3D'
%
%   Output parameters:
%       booelan  - true or false
%
%   TEST_HRTF_EXTRAPOLATION(HRTF_SET) checks if the HRTF exrapolation works
%   correctly.

%*****************************************************************************
% Copyright (c) 2010-2015 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2015 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************

% TODO: add mode to save data as reference data


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);


%% ===== Configuration ===================================================
if strcmp('QU_KEMAR',hrtf_set)
    % QU KEMAR horizontal HRTFs
    disp('Extrapolate QU KEMAR anechoic 3m ...');
    % extrapolation settings
    conf.dimension = '2.5D';
    conf.N = 1024;
    conf.c = 343;
    conf.fs = 44100;
    conf.xref = [0 0 0];
    conf.usetapwin = true;
    conf.tapwinlen = 0.3;
    conf.secondary_sources.geometry = 'circle';
    conf.secondary_sources.center = [0 0 0];
    conf.secondary_sources.size = 3;
    conf.wfs.usehpre = true;
    conf.wfs.hpretype = 'FIR';
    conf.driving_functions = 'default';
    conf.usefracdelay = false;
    conf.fracdelay_method = 'resample';
    conf.ir.useinterpolation = true;
    conf.ir.useoriglength = false;
    conf.showprogress = true;
    % check if HRTF data set is available, download otherwise
    basepath = get_sfs_path();
    hrtf_file = [basepath '/data/HRTFs/QU_KEMAR_anechoic_3m.sofa'];
    if ~exist(hrtf_file,'file')
        url = ['https://raw.githubusercontent.com/sfstoolbox/data/master/', ...
               'HRTFs/QU_KEMAR_anechoic_3m.sofa'];
        download_file(url,hrtf_file);
    end
    % load HRTF data set
    sofa = SOFAload(hrtf_file);
    x0 = SOFAcalculateAPV(sofa);
    % do the extrapolation
    sofa_pw = extrapolate_farfield_hrtfset(sofa,conf);
    % plot the original HRTF data set
    figure;
    imagesc(x0(:,1),1:size(sofa.Data.IR,3),squeeze(sofa.Data.IR(:,1,:))');
    title('QU KEMAR anechoic 3m');
    xlabel('phi / deg');
    % plot the interplated HRTF data set
    figure;
    imagesc(x0(:,1),1:size(sofa_pw.Data.IR,3),squeeze(sofa_pw.Data.IR(:,1,:))');
    title('QU KEMAR anechoic extrapolated');
    xlabel('phi / deg');
    % ILD of both HRTF sets
    ild1 = db(rms(squeeze(sofa.Data.IR(:,1,:))')) - ...
           db(rms(squeeze(sofa.Data.IR(:,2,:))'));
    ild2 = db(rms(squeeze(sofa_pw.Data.IR(:,1,:))')) - ...
           db(rms(squeeze(sofa_pw.Data.IR(:,2,:))'));
    figure;
    plot(x0(:,1),ild1,'-b', ...
         x0(:,1),ild2,'-r');
    legend('original','extrapolated');
    title('Interaural Level Differences');
    xlabel('phi / deg');
    ylabel('Amplitude difference / dB');
elseif strcmp('FABIAN_3D',hrtf_set)
    error(['%s: the FABIAN 3D data set is not publicly available at the ' ...
    'moment.'],upper(mfilename));
    % SEACEN FABIAN 3D HRTFs
    disp('Extrapolate SEACEN FABIAN anechoic ...');
    % extrapolation settings
    conf.dimension = '3D';
    conf.fs = 44100;
    conf.c = 343;
    conf.xref = [0 0 0];
    conf.usetapwin = false;
    conf.tapwinlen = 0.3;
    conf.secondary_sources.geometry = 'sphere';
    conf.N = 1024;
    conf.wfs.usehpre = true;
    conf.wfs.hpretype = 'FIR';
    conf.driving_functions = 'default';
    conf.usefracdelay = false;
    conf.fracdelay_method = 'resample';
    conf.ir.useinterpolation = true;
    conf.ir.useoriglength = false;
    conf.showprogress = true;
    addirspath(conf);
    hrtf_file = 'FABIAN_3d_anechoic.sofa';
    % do the extrapolation
    irs_pw = extrapolate_farfield_hrtfset(hrtf_file,conf);
    %save FABIAN_3D_extrapolated.sofa irs_pw
    % plot the original HRTF data set
    figure;
    imagesc(x0(:,1),1:size(sofa.Data.IR,3),squeeze(sofa.Data.IR(:,1,:))');
    title('SEACEN FABIAN anechoic 1.7m');
    xlabel('phi / deg');
    % plot the interplated HRTF data set
    figure;
    imagesc(x0(:,1),1:size(sofa_pw.Data.IR,3),squeeze(sofa.Data_pw.IR(:,1,:))');
    title('SEACEN FABIAN anechoic extrapolated');
    xlabel('phi / deg');
    % ILD of both HRTF sets
    % ILD of both HRTF sets
    ild1 = db(rms(squeeze(sofa.Data.IR(:,1,:))')) - ...
           db(rms(squeeze(sofa.Data.IR(:,2,:))'));
    ild2 = db(rms(squeeze(sofa_pw.Data.IR(:,1,:))')) - ...
           db(rms(squeeze(sofa_pw.Data.IR(:,2,:))'));
    figure;
    plot(x0(:,1),ild1,'-b', ...
         x0(:,1),ild2,'-r');
    legend('original','extrapolated');
    title('Interaural Level Differences');
    xlabel('phi / deg');
    ylabel('Amplitude difference / dB');
end
