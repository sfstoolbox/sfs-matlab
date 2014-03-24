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
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
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
    hrtf_file = [basepath '/data/HRTFs/QU_KEMAR_anechoic_3m.mat'];
    if ~exist(hrtf_file,'file')
        url = ['https://dev.qu.tu-berlin.de/projects/measurements/repository/', ...
            'raw/2010-11-kemar-anechoic/mat/QU_KEMAR_anechoic_3m.mat'];
        download_file(url,hrtf_file);
    end
    % load HRTF data set
    irs = read_irs(hrtf_file,conf);
    % do the extrapolation
    irs_pw = extrapolate_farfield_hrtfset(irs,conf);
    % plot the original HRTF data set
    figure;
    imagesc(deg(irs.apparent_azimuth),1:size(irs.left,1),irs.left);
    title('QU KEMAR anechoic 3m');
    xlabel('phi / deg');
    % plot the interplated HRTF data set
    figure;
    imagesc(deg(irs_pw.apparent_azimuth),1:size(irs_pw.left,1),irs_pw.left);
    title('QU KEMAR anechoic extrapolated');
    xlabel('phi / deg');
    % ILD of both HRTF sets
    ild1 = interaural_level_difference(irs.left,irs.right);
    ild2 = interaural_level_difference(irs_pw.left,irs_pw.right);
    figure;
    plot(deg(irs.apparent_azimuth),ild1,'-b', ...
         deg(irs.apparent_azimuth),ild2,'-r');
    legend('original','extrapolated');
    title('Interaural Level Differences');
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
    hrtf_file = 'FABIAN_3d_anechoic.mat';
    % load HRTF data set
    irs = read_irs(hrtf_file,conf);
    % do the extrapolation
    irs_pw = extrapolate_farfield_hrtfset(irs,conf);
    %save FABIAN_3D_extrapolated.mat irs_pw
    % plot the original HRTF data set
    figure;
    imagesc(deg(irs.apparent_azimuth),1:size(irs.left,1),irs.left);
    title('SEACEN FABIAN anechoic 1.7m');
    xlabel('phi / deg');
    % plot the interplated HRTF data set
    figure;
    imagesc(deg(irs_pw.apparent_azimuth),1:size(irs_pw.left,1),irs_pw.left);
    title('SEACEN FABIAN anechoic extrapolated');
    xlabel('phi / deg');
    % ILD of both HRTF sets
    ild1 = interaural_level_difference(irs.left,irs.right);
    ild2 = interaural_level_difference(irs_pw.left,irs_pw.right);
    figure;
    plot(deg(irs.apparent_azimuth),ild1,'-b', ...
         deg(irs.apparent_azimuth),ild2,'-r');
    legend('original','extrapolated');
    title('Interaural Level Differences');
end
