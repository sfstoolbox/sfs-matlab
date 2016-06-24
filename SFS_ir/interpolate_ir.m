function [ir_new,x0_new] = interpolate_ir(ir,x0,xs,conf)
%INTERPOLATE_IR interpolates three given IRs for the given angle
%
%   Usage: [ir,x0] = interpolate_ir(ir,x0,xs,conf)
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
%       ir      - impulse response for the given position [1 C N]
%       x0      - position corresponding to the returned impulse response
%
%   INTERPOLATE_IR(ir,x0,xs,conf) interpolates the two to three given impulse
%   responses from ir with their corresponding angles x0 for the given angles
%   xs and returns an interpolated impulse response.
%   Note, that the given parameter are not checked if they have all the correct
%   dimensions in order to save computational time, because this function could
%   be called quite often.
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
if useinterpolation && ~isfield(conf.ir,'interpolationmethod')
    warning('SFS:irs_intpolmethod',...
        'no interpolation method provided, will use method ''simple''.');
    interpolationmethod = 'simple';
else
    interpolationmethod = conf.ir.interpolationmethod;
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
                TF = fft(ir,[],3);
                % Calculate interpolation for magnitude and phase separately
                mag = abs(TF);
                pha = unwrap(angle(TF),[],3);
                % Calculate interpolation only for the first half of the spectrum
                idx_half = floor(size(mag,3)/2)+1;
                mag_new(1,1,:) = interpolation(squeeze(mag(1:2,1,1:idx_half))',...
                    x0(:,1:2),xs);
                mag_new(1,2,:) = interpolation(squeeze(mag(1:2,2,1:idx_half))',...
                    x0(:,1:2),xs);
                pha_new(1,1,:) = interpolation(squeeze(pha(1:2,1,1:idx_half))',...
                    x0(:,1:2),xs);
                pha_new(1,2,:) = interpolation(squeeze(pha(1:2,2,1:idx_half))',...
                    x0(:,1:2),xs);
                ir_new(1,1,:) = ifft(squeeze(mag_new(1,1,:)...
                    .*exp(1i*pha_new(1,1,:))),size(mag,3),'symmetric');
                ir_new(1,2,:) = ifft(squeeze(mag_new(1,2,:)...
                    .*exp(1i*pha_new(1,1,:))),size(mag,3),'symmetric');
            otherwise 
                error('%s: %s is an unknown interpolation method.',...
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
