function irs_pw = extrapolate_farfield_hrtfset(irs,conf)
%EXTRAPOLATE_FARFIELD_HRTFSET far-field extrapolation of a given HRTF dataset
%
%   Usage: irs = extrapolate_farfield_hrtfset(irs,[conf])
%
%   Input parameters:
%       irs     - IR data set for the virtual secondary sources
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       irs     - IR data set extra polated to conation plane wave IRs
%
%   EXTRAPOLATE_FARFIELD_HRTFSET(IRS) generates a far-field extrapolated set of
%   impulse responses, using the given irs set. Far-field means that the
%   resulting impulse responses are plane waves. The extrapolation is done via
%   WFS.
%
%   References:
%       S. Spors and J. Ahrens (2011) - "Generation of far-field head-related
%       transfer functions using sound field synthesis", In German Annual
%       Conference on Acoustics (DAGA).
%
%   see also: get_ir, driving_function_imp_wfs

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


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
check_irs(irs);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ===================================================
fs = conf.fs;
dimension = conf.dimension;
N = size(irs.left,1);
showprogress = conf.showprogress;
conf.ir.usehcomp = false;


%% ===== Variables ======================================================
nls = length(irs.apparent_azimuth);
if length(irs.distance)==1
    R = repmat(irs.distance,nls,1);
else
    R = irs.distance';
end
phi = irs.apparent_azimuth';
theta = irs.apparent_elevation';
conf.secondary_sources.number = nls;
conf.secondary_sources.x0 = zeros(nls,7);
[conf.secondary_sources.x0(:,1), ...
 conf.secondary_sources.x0(:,2), ...
 conf.secondary_sources.x0(:,3)] = sph2cart(phi,theta,R);
conf.secondary_sources.x0(:,4:6) = ...
    direction_vector(conf.secondary_sources.x0(:,1:3),repmat(conf.xref,nls,1));
% weights
if strcmp('3D',dimension)
    % use rectangular grid to get a first approximation of the grid weights
    % R^2 * cos(theta) is the integrational weight for integration on a sphere
    conf.secondary_sources.x0(:,7) = ...
        weights_for_points_on_a_sphere_rectangle(phi,theta) .* ...
        R.^2 .* cos(theta);
else
    conf.secondary_sources.x0(:,7) = ones(nls,1);
end
% check if we have a 2D or 3D secondary source setup
if any(irs.apparent_elevation-irs.apparent_elevation(1)) && ...
    any(irs.apparent_azimuth-irs.apparent_azimuth(1))
    % 3D case
    conf.usetapwin = false;
    if ~strcmp('3D',dimension)
        warning(['You have a 3D HRTF data set, but are not using ', ...
            'conf.dimension="3D".']);
    end
else
    if strcmp('3D',dimension)
        warning(['You are using a 2D HRTF data set, but are using ', ...
            'conf.dimension="3D".']);
    end
end
if strcmp('2.5D',dimension)
    % Apply a amplitude correction, due to 2.5D. This will result in a correct
    % reproduced ILD in the resulting impulse responses (see, Spors 2011)
    amplitude_correction = -1.7 * sin(irs.apparent_azimuth);
else
    amplitude_correction = zeros(size(irs.apparent_azimuth));
end


%% ===== Computation =====================================================
% get virtual secondary source positions
x0_all = secondary_source_positions(conf);
conf.wfs.hpreflow = 50;
conf.wfs.hprefhigh = aliasing_frequency(x0_all,conf);

% Initialize new irs set
irs_pw = irs;
irs_pw.description = 'Extrapolated HRTF set containing plane waves';
irs_pw.left = zeros(N,nls);
irs_pw.right = zeros(N,nls);
irs_pw.distance = Inf;


% Get all HRTFs for the secondary source positions
ir_all = zeros(size(irs.left,1),2,nls);
for ii=1:nls
    ir_all(:,:,ii) = get_ir(irs,x0_all(ii,1:3),'cartesian',conf);
end

% Generate a irs set for all given angles
for ii=1:nls

    % show progress
    if showprogress, progress_bar(ii,nls); end;

    % direction of plane wave
    [xs(1),xs(2),xs(3)] = sph2cart(phi(ii),theta(ii),R(ii));
    xs = -xs;

    % calculate active virtual speakers
    [x0,idx] = secondary_source_selection(x0_all,xs,'pw');
    ir = ir_all(:,:,idx);
    % apply tapering window
    x0 = secondary_source_tapering(x0,conf);

    % get driving signals, temporarely deactivate WFS pre-filter, because it
    % will be applied once at the end
    tmp_usehpre = conf.wfs.usehpre;
    conf.wfs.usehpre = false;
    [~,delay,weight] = driving_function_imp_wfs(x0,xs,'pw',conf);
    conf.wfs.usehpre = tmp_usehpre;
    % delay in samples
    delay = delay.*fs;
    % sum up contributions from individual virtual speakers
    irs_pw.left(:,ii) = sum(delayline(squeeze(ir(:,1,:)),delay,weight,conf),2);
    irs_pw.right(:,ii) = sum(delayline(squeeze(ir(:,2,:)),delay,weight,conf),2);
    irs_pw.left(:,ii) = irs_pw.left(:,ii)/10^(amplitude_correction(ii)/20);
    irs_pw.right(:,ii) = irs_pw.right(:,ii)/10^(-amplitude_correction(ii)/20);

end

%% ===== Pre-equalization ===============================================
irs_pw.left = wfs_preequalization(irs_pw.left,conf);
irs_pw.right = wfs_preequalization(irs_pw.right,conf);
