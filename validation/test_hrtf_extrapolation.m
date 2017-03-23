function status = test_hrtf_extrapolation(modus,hrtf_set)
%TEST_HRTF_EXTRAPOLATION tests the HRTF extrapolation functions
%
%   Usage: status = test_hrtf_extrapolation(modus,hrtf_set)
%
%   Input parameters:
%       modus    - 0: numerical
%                  1: visual
%       hrtf_set - HRTFs for testing extrapolation. Can be
%                    'QU_KEMAR'
%                    'FABIAN_3D'
%
%   Output parameters:
%       status - true or false
%
%   TEST_HRTF_EXTRAPOLATION(modus,hrtf_set) checks if the HRTF exrapolation works
%   correctly.

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


% TODO: add mode to save data as reference data
status = false;


%% ===== Checking of input  parameters ===================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);


%% ===== Configuration ===================================================
conf = SFS_config;
conf.showprogress = true;
if strcmp('QU_KEMAR',hrtf_set)
    % QU KEMAR horizontal HRTFs
    disp('Extrapolate QU KEMAR anechoic 3m ...');
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
    if modus
        figure;
        imagesc(x0(:,1),1:size(sofa.Data.IR,3),squeeze(sofa.Data.IR(:,1,:))');
        title('QU KEMAR anechoic 3m');
        xlabel('phi / deg');
        % plot the interplated HRTF data set
        figure;
        imagesc(x0(:,1),1:size(sofa_pw.Data.IR,3),squeeze(sofa_pw.Data.IR(:,1,:))');
        title('QU KEMAR anechoic extrapolated');
        xlabel('phi / deg');
    end
    % ILD of both HRTF sets
    ild1 = db(rms(squeeze(sofa.Data.IR(:,1,:))')) - ...
           db(rms(squeeze(sofa.Data.IR(:,2,:))'));
    ild2 = db(rms(squeeze(sofa_pw.Data.IR(:,1,:))')) - ...
           db(rms(squeeze(sofa_pw.Data.IR(:,2,:))'));
    if modus
        figure;
        plot(x0(:,1),ild1,'-b', ...
             x0(:,1),ild2,'-r');
        legend('original','extrapolated');
        title('Interaural Level Differences');
        xlabel('phi / deg');
        ylabel('Amplitude difference / dB');
    end
elseif strcmp('FABIAN_3D',hrtf_set)
    error(['%s: the FABIAN 3D data set is not publicly available at the ' ...
    'moment.'],upper(mfilename));
    % SEACEN FABIAN 3D HRTFs
    disp('Extrapolate SEACEN FABIAN anechoic ...');
    % extrapolation settings
    conf.dimension = '3D';
    conf.usetapwin = false;
    conf.secondary_sources.geometry = 'sphere';
    hrtf_file = 'FABIAN_3d_anechoic.sofa';
    % do the extrapolation
    irs_pw = extrapolate_farfield_hrtfset(hrtf_file,conf);
    %save FABIAN_3D_extrapolated.sofa irs_pw
    % plot the original HRTF data set
    if modus
        figure;
        imagesc(x0(:,1),1:size(sofa.Data.IR,3),squeeze(sofa.Data.IR(:,1,:))');
        title('SEACEN FABIAN anechoic 1.7m');
        xlabel('phi / deg');
        % plot the interplated HRTF data set
        figure;
        imagesc(x0(:,1),1:size(sofa_pw.Data.IR,3),squeeze(sofa.Data_pw.IR(:,1,:))');
        title('SEACEN FABIAN anechoic extrapolated');
        xlabel('phi / deg');
    end
    % ILD of both HRTF sets
    % ILD of both HRTF sets
    ild1 = db(rms(squeeze(sofa.Data.IR(:,1,:))')) - ...
           db(rms(squeeze(sofa.Data.IR(:,2,:))'));
    ild2 = db(rms(squeeze(sofa_pw.Data.IR(:,1,:))')) - ...
           db(rms(squeeze(sofa_pw.Data.IR(:,2,:))'));
    if modus
        figure;
        plot(x0(:,1),ild1,'-b', ...
             x0(:,1),ild2,'-r');
        legend('original','extrapolated');
        title('Interaural Level Differences');
        xlabel('phi / deg');
        ylabel('Amplitude difference / dB');
    end
end


status = true;
