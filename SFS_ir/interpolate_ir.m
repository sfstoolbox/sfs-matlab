function [ir_new,x0_new] = interpolate_ir(ir,x0,xs,conf)
%INTERPOLATE_IR interpolates three given IRs for the given angle
%
%   Usage: [ir_new,x0_new] = interpolate_ir(ir,x0,xs,conf)
%
%   Input parameters:
%       ir      - matrix containing impulse responses in the form [M C N], where
%                     M ... Number of measurements (2<=M<=3)
%                     C ... Number of channels
%                     N ... Number of samples
%       x0      - matrix containing positions of single impulse
%                 responses [2 M] / (rad, rad)
%       xs      - desired position after interpolation [2 1] / (rad, rad)
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       ir_new  - impulse response for the given position [1 C N]
%       x0_new  - position corresponding to the returned impulse response
%
%   INTERPOLATE_IR(ir,x0,xs,conf) interpolates the two to three given impulse
%   responses from ir with their corresponding angles x0 for the given angles
%   xs and returns an interpolated impulse response.
%	For the 1D case, the interpolation method differs depending on the setting
%	of conf.ir.interpolationmethod:
%     'simple'      - Interpolation in the time domain performed samplewise.
%	                  This does not heed the times of arrival of the impulse
%                     responses.
%     'freqdomain'  - Interpolation in the frequency domain performed separately
%                     for magnitude and phase.
%   Note, that the given parameters are not checked if they have all the correct
%   dimensions in order to save computational time, because this function could
%   be called quite often.
%
%	References:
%		K. Hartung, J. Braasch, S. J. Sterbing (1999) - "Comparison of different
%		methods for the interpolation of head-related transfer functions".
%		Proc. of the 16th AES Conf.
%		K. Itoh (1982) - "Analysis of the phase unwrapping algorithm". Applied
%		Optics 21(14), 2470
%
%   See also: get_ir, interpolation

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Developers                             *
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


%% ===== Checking of input parameters ===================================
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);


%% ===== Configuration ==================================================
useinterpolation = conf.ir.useinterpolation;
% Check for old configuration
if useinterpolation 
    if ~isfield(conf.ir,'interpolationmethod')
        warning('SFS:irs_intpolmethod',...
            'no interpolation method provided, will use method ''simple''.');
        interpolationmethod = 'simple';
    else
        interpolationmethod = conf.ir.interpolationmethod;
    end
end
% Precision of the wanted angle. If an impulse response within the given
% precision could be found no interpolation is applied.
prec = 0.001; % ~ 0.05 deg


%% ===== Computation ====================================================
ir_new = ir(1,:,:);
x0_new = xs;
% Check if we have found directly the desired point or have to interpolate
% between different impulse responses
if norm(x0(:,1)-xs)<prec || ~useinterpolation || size(x0,2)==1
    % Return the first nearest neighbour
    x0_new = x0(:,1);
    return;
else
    % === IR interpolation ===
    % Check if we have to interpolate in one or two dimensions
    if norm(x0(1,1)-x0(1,2))<prec || norm(x0(2,1)-x0(2,2))<prec
        % --- 1D interpolation ---
        warning('SFS:irs_intpol',...
            ['doing 1D IR interpolation between (%.1f,%.1f) deg ', ...
             'and (%.1f,%.1f) deg.'], ...
            deg(x0(1,1)), deg(x0(2,1)), ...
            deg(x0(1,2)), deg(x0(2,2)));
        switch interpolationmethod
        case 'simple'
            ir_new(1,1,:) = interpolation(squeeze(ir(1:2,1,:))',x0(:,1:2),xs);
            ir_new(1,2,:) = interpolation(squeeze(ir(1:2,2,:))',x0(:,1:2),xs);
        case 'freqdomain'
            % see Itoh (1982), Hartung et al. (1999)
            %
            % Upsample to avoid phase aliasing in unwrapping of phase
            TF = fft(ir,4*size(ir,3),3);
            % Magnitude and phase will be interpolated separately
            magnitude = abs(TF);
            phase = unwrap(angle(TF),[],3);
            % Calculate interpolation only for the first half of the spectrum
            % and only for original bins
            idx_half = floor(size(TF,3)/2)+1;
            magnitude_new(1,:) = interpolation(...
                squeeze(magnitude(1:2,1,1:4:idx_half))',x0(:,1:2),xs);
            magnitude_new(2,:) = interpolation(...
                squeeze(magnitude(1:2,2,1:4:idx_half))',x0(:,1:2),xs);
            phase_new(1,:) = interpolation(...
                squeeze(phase(1:2,1,1:4:idx_half))',x0(:,1:2),xs);
            phase_new(2,:) = interpolation(...
                squeeze(phase(1:2,2,1:4:idx_half))',x0(:,1:2),xs);
            % Calculate interpolated impulse response from magnitude and phase
            ir_new(1,1,:) = ifft(magnitude_new(1,:) ...
                .* exp(1i*phase_new(1,:)),size(ir,3),'symmetric');
            ir_new(1,2,:) = ifft(magnitude_new(2,:)...
                .* exp(1i*phase_new(2,:)),size(ir,3),'symmetric');
        otherwise
            error('%s: %s is an unknown interpolation method.', ...
                upper(mfilename),interpolationmethod);
        end
    else
        % --- 2D interpolation ---
        warning('SFS:irs_intpol3D',...
            ['doing 2D IR interpolation between (%.1f,%.1f) deg, ', ...
             '(%.1f,%.1f) deg and (%.1f,%.1f) deg.'], ...
            deg(x0(1,1)), deg(x0(2,1)), ...
            deg(x0(1,2)), deg(x0(2,2)), ...
            deg(x0(1,3)), deg(x0(2,3)));
        ir_new(1,1,:) = interpolation(squeeze(ir(1:3,1,:))',x0(:,1:3),xs);
        ir_new(1,2,:) = interpolation(squeeze(ir(1:3,2,:))',x0(:,1:3),xs);
    end
end
