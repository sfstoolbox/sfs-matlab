function sofa_pw = extrapolate_farfield_hrtfset(sofa,conf)
%EXTRAPOLATE_FARFIELD_HRTFSET far-field extrapolation of a given HRTF dataset
%
%   Usage: sofa = extrapolate_farfield_hrtfset(sofa,conf)
%
%   Input parameters:
%       sofa    - IR data set for the virtual secondary sources
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       sofa    - IR data set extra polated to conation plane wave IRs
%
%   EXTRAPOLATE_FARFIELD_HRTFSET(SOFA,conf) generates a far-field extrapolated
%   set of impulse responses, using the given irs set. Far-field means that the
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


%% ===== Checking of input  parameters ==================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargstruct(sofa,conf);


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
    direction_vector(conf.secondary_sources.x0,repmat(conf.xref,nls,1));
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
warning('off','SFS:get_ir'); % Disable warning for short N
for ii=1:nls
    ir_all(ii,:,:) = get_ir(sofa,X,head_orientation, ...
                            x0_all(ii,1:3),'cartesian',conf)';
end
warning('on','SFS:get_ir');

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
