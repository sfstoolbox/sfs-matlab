function ir = interpolate_ir(ir,weights,conf)
%INTERPOLATE_IR interpolates the given impulse responses according to their weights
%
%   Usage: ir = interpolate_ir(ir,weights,conf)
%
%   Input parameters:
%       ir           - matrix containing impulse responses in the form [M C N], where
%                          M ... Number of measurements
%                          C ... Number of channels
%                          N ... Number of samples
%       weights      - M weights for impulse reponses [M 1]
%       conf         - configuration struct (see SFS_config)
%
%   Output parameters:
%       ir           - impulse response for the given position [1 C N]
%
%   INTERPOLATE_IR(ir,weights,conf) interpolates the given impulse responses
%   by applying the given weights and returns the interpolated impulse response.
%   Only impulse responses with weights larger that the precision prec=0.001 will
%   be used.
%	The interpolation method differs depending on the setting of
%	conf.ir.interpolationmethod:
%     'simple'      - Interpolation in the time domain performed samplewise.
%	                  This does not heed the times of arrival of the impulse
%                     responses.
%     'freqdomain'  - Interpolation in the frequency domain performed separately
%                     for magnitude and phase.
%     'timedomain'  - Interpolation in the time domain with cross-correlation for
%                     estimation of time of arrival (TOA) differences
%   Note that the given parameters are not checked if they all have the correct
%   dimensions in order to save computational time, because this function could
%   be called quite often.
%
%   See also: get_ir, interpolation
%
%	References:
%		Hartung, Braasch, Sterbing (1999) - "Comparison of different methods for
%		the interpolation of head-related transfer functions", 16th Conference
%		of the Audio Engineering Society, Paper 16-028,
%		http://www.aes.org/e-lib/browse.cfm?elib=8026
%
%		Itoh (1982) - "Analysis of the phase unwrapping algorithm", Applied
%		Optics 21(14), p. 2470, https://doi.org/10.1364/AO.21.002470

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2019 SFS Toolbox Developers                             *
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
% https://sfs.readthedocs.io                            sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input parameters ===================================
nargmin = 3;
nargmax = 3;
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
weights = weights(weights>=prec);

% === IR interpolation ===
if useinterpolation && length(weights)>1
    switch interpolationmethod
    case 'simple'
        ir = sum(bsxfun(@times,ir,weights),1);
    case 'freqdomain'
        % See Itoh (1982), Hartung et al. (1999)
        %
        % Upsample to avoid phase aliasing in unwrapping of phase
        TF = fft(ir,4*size(ir,3),3);
        % Magnitude and phase will be interpolated separately
        idx_half = 2*size(ir,3)+1;  % index for first half of spectrum
        magnitude = abs(TF(:,:,1:idx_half));
        phase = unwrap(angle(TF(:,:,1:idx_half)),[],3);
        % Mirror magnitude and phase
        magnitude = cat(3,magnitude,magnitude(:,:,idx_half-1:-1:2));
        phase = cat(3,phase,-phase(:,:,idx_half-1:-1:2));
        % Calculate interpolation of the spectrum and downsample
        magnitude = sum(bsxfun(@times,magnitude(:,:,1:4:end),weights),1);
        phase = sum(bsxfun(@times,phase(:,:,1:4:end),weights),1);
        % Calculate interpolated impulse response from new magnitude and phase
        ir = ifft(magnitude.*exp(1i*phase),[],3);
        % Avoid round-off errors
        ir = real(ir);
    case 'timedomain'
        % See Hartung et al. (1999)
        %
        % Loop over channels
        ir_results = zeros(1,size(ir,2),size(ir,3));
        for n = 1:size(ir,2)
            % Determine TOA differences
            ir_upsampled = resample(squeeze(ir(:,n,:))',10,1);
            for nn = 1:size(ir,1)
                [~,idx_max(nn)] = ...
                    max(xcorr(ir_upsampled(:,1),ir_upsampled(:,nn)));
            end
            N = size(ir_upsampled,1);
            k = -(N-1):(N-1);
            shifts = k(idx_max);
            
            % Shift all irs back to align with last ir
            max_shift = max(shifts)-min(shifts);
            TOA_diff_to_last = shifts-min(shifts);
            ir_shifted = zeros(N+max_shift,size(ir,1));
            for nn = 1:size(ir,1)
                ir_shifted(:,nn) = [zeros(TOA_diff_to_last(nn),1); ...
                    ir_upsampled(:,nn); ...
                    zeros(max_shift-TOA_diff_to_last(nn),1)];
            end
            
            % Interpolate aligned irs
            ir_interp = sum(bsxfun(@times,ir_shifted,weights'),2);
            
            % Interpolate TOA differences and shift interpolated ir accordingly
            TOA_diff_interp = round(TOA_diff_to_last*weights);
            ir_interp = ...
                ir_interp(TOA_diff_interp+1:end-(max_shift-TOA_diff_interp));

            % Downsample and collect results per channel
            ir_results(1,n,:) = resample(ir_interp,1,10);
        end
        ir = ir_results;
    otherwise
        error('%s: %s is an unknown interpolation method.', ...
            upper(mfilename),interpolationmethod);
    end
end
