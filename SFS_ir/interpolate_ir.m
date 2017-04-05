function [ir_new,weights_new,x0_new] = interpolate_ir(ir,weights,x0,conf)
%INTERPOLATE_IR interpolates the given impulse responses according to their weights
%
%   Usage: [ir_new,weights_new,x0_new] = interpolate_ir(ir,weights,x0,conf)
%
%   Input parameters:
%       ir           - matrix containing impulse responses in the form [M C N], where
%                          M ... Number of measurements
%                          C ... Number of channels
%                          N ... Number of samples
%       weights      - M weights for impulse reponses
%       x0           - M positions corresponding to given impulse responses
%       conf         - configuration struct (see SFS_config)
%
%   Output parameters:
%       ir_new       - impulse response for the given position [1 C N]
%       weights_new  - weights corresponding to impulse responses used for
%                      interpolation
%       x0_new       - position corresponding to impulse responses used for
%                      interpolation
%
%   INTERPOLATE_IR(ir,x0,xs,conf) interpolates the given impulse responses by
%   applying the given weights and returns the interpolated impulse response as
%   well as its corresponding position. Only impulse responses with weights larger
%   that the precision prec=0.001 will be used.
%	The interpolation method differs depending on the setting of
%	conf.ir.interpolationmethod:
%     'simple'      - Interpolation in the time domain performed samplewise.
%	                  This does not heed the times of arrival of the impulse
%                     responses.
%     'freqdomain'  - Interpolation in the frequency domain performed separately
%                     for magnitude and phase.
%   Note that the given parameters are not checked if they all have the correct
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

% Precision of the weights. Impulse responses with smaller weights are left out.
prec = 0.001;


%% ===== Computation ====================================================
% Leave out impulse responses with weights smaller than prec
ir = ir(weights>=prec,:,:);
weights_new = weights(weights>=prec);
x0_new = x0(weights>=prec,:);

% === IR interpolation ===
if ~useinterpolation || length(weights_new)==1
    ir_new = ir;
elseif useinterpolation
    switch interpolationmethod
    case 'simple'
        ir_new = sum(bsxfun(@times,ir,weights),1);
    case 'freqdomain'
        % See Itoh (1982), Hartung et al. (1999)
        %
        % Upsample to avoid phase aliasing in unwrapping of phase
        TF = fft(ir,4*size(ir,3),3);
        % Magnitude and phase will be interpolated separately
        magnitude = abs(TF);
        phase = unwrap(angle(TF),[],3);
        % Calculate interpolation only for the first half of the spectrum
        % and only for original bins
        idx_half = floor(size(TF,3)/2)+1;
        magnitude_new = sum(bsxfun(@times,magnitude(:,:,1:4:idx_half),weights),1);
        phase_new = sum(bsxfun(@times,phase(:,:,1:4:idx_half),weights),1);
        % Calculate interpolated impulse response from new magnitude and phase
        ir_new = ifft(magnitude_new.*exp(1i*phase_new),size(ir,3),3,'symmetric');
    otherwise
        error('%s: %s is an unknown interpolation method.', ...
            upper(mfilename),interpolationmethod);
    end   
end
