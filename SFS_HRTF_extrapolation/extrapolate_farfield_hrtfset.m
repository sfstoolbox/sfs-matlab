function sofa_pw = extrapolate_farfield_hrtfset(sofa,conf)
%EXTRAPOLATE_FARFIELD_HRTFSET far-field extrapolation of a given HRTF dataset
%
%   Usage: sofa = extrapolate_farfield_hrtfset(sofa,[conf])
%
%   Input parameters:
%       sofa    - IR data set for the virtual secondary sources
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       sofa    - IR data set extra polated to conation plane wave IRs
%
%   EXTRAPOLATE_FARFIELD_HRTFSET(SOFA) generates a far-field extrapolated set of
%   impulse responses, using the given irs set. Far-field means that the
%   resulting impulse responses are plane waves. The extrapolation is done via
%   WFS.
%
%   References:
%       S. Spors and J. Ahrens (2011) - "Generation of far-field head-related
%       transfer functions using sound field synthesis", In German Annual
%       Conference on Acoustics (DAGA).
%
%   See also: get_ir, driving_function_imp_wfs

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


%% ===== Checking of input  parameters ==================================
nargmin = 1;
nargmax = 2;
narginchk(nargmin,nargmax);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end
isargstruct(sofa);


%% ===== Configuration ===================================================
fs = conf.fs;
dimension = conf.dimension;
showprogress = conf.showprogress;
conf.ir.usehcomp = false;


%% ===== Variables ======================================================
[nls,~,N] = size(sofa.Data.IR);
APV = SOFAcalculateAPV(sofa);
phi = rad(APV(:,1));
theta = rad(APV(:,2));
R = APV(:,3);
conf.secondary_sources.number = nls;
[conf.secondary_sources.x0(:,1), ...
 conf.secondary_sources.x0(:,2), ...
 conf.secondary_sources.x0(:,3)] = sph2cart(phi,theta,R);
conf.secondary_sources.x0(:,4:6) = ...
    direction_vector(conf.secondary_sources.x0, ...
                     repmat(conf.secondary_sources.center,nls,1));
% Weights
if strcmp('3D',dimension)
    % Use rectangular grid to get a first approximation of the grid weights
    % R^2 * cos(theta) is the integrational weight for integration on a sphere
    conf.secondary_sources.x0(:,7) = ...
        weights_for_points_on_a_sphere_rectangle(phi,theta) .* ...
        R.^2 .* cos(theta);
else
    conf.secondary_sources.x0(:,7) = ones(nls,1);
end
% Check if we have a 2D or 3D secondary source setup
if any(phi-phi(1)>eps('single')) && any(theta-theta(1)>eps('single'))
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
    amplitude_correction = -1.7 * sin(phi);
else
    amplitude_correction = zeros(nls,1);
end


%% ===== Computation =====================================================
% Get virtual secondary source positions
x0_all = secondary_source_positions(conf);
conf.wfs.hpreflow = 50;
conf.wfs.hprefhigh = aliasing_frequency(x0_all,conf);

% Initialize new irs set
sofa_pw = sofa;
sofa_pw.GLOBAL_Comment = 'Extrapolated HRTF set containing plane waves';
sofa_pw.Data.IR = zeros(nls,2,N);
sofa_pw.SourcePosition = [Inf 0 0];

% Get all HRTFs for the secondary source positions
ir_all = zeros(nls,2,N);
X = SOFAconvertCoordinates(sofa.ListenerPosition, ...
                           sofa.ListenerPosition_Type, ...
                           'cartesian');
head_orientation = SOFAconvertCoordinates(sofa.ListenerView, ...
                                          sofa.ListenerView_Type, ...
                                          'spherical');
head_orientation = rad(head_orientation(1,1:2));
for ii=1:nls
    ir_all(ii,:,:) = get_ir(sofa,X,head_orientation, ...
                            x0_all(ii,1:3),'cartesian',conf)';
end

% Generate a impulse response set for all given angles
for ii=1:nls

    % Show progress
    if showprogress, progress_bar(ii,nls); end;

    % Direction of plane wave
    [xs(1),xs(2),xs(3)] = sph2cart(phi(ii),theta(ii),R(ii));
    xs = -xs;

    % Calculate active virtual speakers
    [x0,idx] = secondary_source_selection(x0_all,xs,'pw');
    ir = ir_all(idx,:,:);
    % Apply tapering window
    x0 = secondary_source_tapering(x0,conf);

    % Get driving signals, temporarely deactivate WFS pre-filter, because it
    % will be applied once at the end
    tmp_usehpre = conf.wfs.usehpre;
    conf.wfs.usehpre = false;
    [~,delay,weight] = driving_function_imp_wfs(x0,xs,'pw',conf);
    conf.wfs.usehpre = tmp_usehpre;
    % Delay in samples
    delay = delay.*fs;
    % Sum up contributions from individual virtual speakers
    for jj=1:size(x0,1)
        % Delay and weight HRTFs
        sofa_pw.Data.IR(ii,:,:) = squeeze(sofa_pw.Data.IR(ii,:,:)) + ...
            delayline(squeeze(ir(jj,:,:))',delay(jj),weight(jj),conf)';
    end
    sofa_pw.Data.IR(ii,1,:) = sofa_pw.Data.IR(ii,1,:)/10^(amplitude_correction(ii)/20);
    sofa_pw.Data.IR(ii,2,:) = sofa_pw.Data.IR(ii,2,:)/10^(-amplitude_correction(ii)/20);

end

%% ===== Pre-equalization ===============================================
sofa_pw.Data.IR(:,1,:) = wfs_preequalization(squeeze(sofa_pw.Data.IR(:,1,:))',conf)';
sofa_pw.Data.IR(:,2,:) = wfs_preequalization(squeeze(sofa_pw.Data.IR(:,2,:))',conf)';
