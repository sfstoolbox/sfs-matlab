function [win,varargout] = modal_weighting(order,ndtft,conf)
%MODAL_WEIGHTING computes weighting window for modal coefficients
%
%   Usage: [win,Win,Phi] = modal_weighting(order,[ndtft],conf)
%
%   Input parameters:
%       order       - half width of weighting window / 1
%       ndtft       - number of bins for inverse discrete-time Fourier transform
%                     (DTFT) / 1 (optional, default: 2*order+1)
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       win         - the window w_n in the discrete domain, only positive n
%                     (length = order+1)
%       Win         - the inverse DTFT of w_n (length = ndtft)
%       Phi         - corresponding angle the inverse DTFT of w_n
%
%   MODAL_WEIGHTING(order,ndtft,conf) calculates a weighting window for the
%   modal band limitation applied in NFC-HOA and LSFS-SBL. The window type is
%   configured in conf.modal_window. Its default setting is a simple
%   rectangular window, for other options have a look into SFS_config. The
%   windows may be different for the 2D/2.5D and the 3D case.
%
%   References:
%   	Kaiser, J., & Schafer, R. (1980) - "On the use of the I0-sinh window
%           for spectrum analysis", IEEE Transactions on Acoustics, Speech, and
%           Signal Processing
%     Daniel, J., Rault, J.-B., Polack, J.-D. (1998) "Ambisonics Encoding of 
%           Other Audio Formats for Multiple Listening Conditions", Proc. of 
%           105th Aud. Eng. Soc. Conv.
%     Zotter, F. & Frank, M. (2012) - "All-Round Ambisonic Panning and
%           Decoding", J. Aud. Eng. Soc. 
%     Van Trees, H. L. (2004) - "Optimum Array Processing", John Wiley & Sons.
%
%   See also: driving_function_imp_nfchoa, driving_function_mono_nfchoa

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


%% ===== Checking input parameters =======================================
nargmin = 2;
nargmax = 3;
narginchk(nargmin,nargmax);
isargpositivescalar(order);
if nargin<nargmax
    conf = ndtft;
    ndtft = 2*order + 1;
end
isargpositivescalar(ndtft);
isargstruct(conf);


%% ===== Configuration ===================================================
wtype = conf.modal_window;


%% ===== Computation =====================================================
switch wtype
case 'rect'
    % === Rectangular Window =========================================
    win = ones(1,order+1);
case 'max-rE'
    % === max-rE window ==============================================
    if any( strcmp(dimension, {'2D', '2.5D'}) )
        % The two-dimensional max-rE window is basically a modified cosine 
        % window, which yields zero for m=order+1 instead of m=order. Hence its
        % last value is not zero. See Daniel (1998), Eq. (44)
        win = cos(pi./2.*(0:order)/(order+1));
    else
        % Approximate solution for the three-dimensional max-rE optimisation 
        % problem. See Zotter (2012), Eq. (10)
        win = zeros(1,order+1);
        for n=0:order
            win(n+1) = asslegendre(n,0,cosd(137.9/(order+1.51)));
        end         
    end    
case {'kaiser', 'kaiser-bessel'}
    % === Kaiser-Bessel window =======================================
    % Approximation of the slepian window using modified bessel
    % function of zeroth order, see Kaiser (1980)
    beta = conf.modal_window_parameter * pi;
    win = besseli(0,beta*sqrt(1-((0:order)./order).^2)) ./ ...
          besseli(0,beta);
case 'tukey'
    % === modified Tukey window ======================================
    % The original tukey is sometimes referred to as the tapered cosine window.
    % It yields unity for m <= alpha*order. For m > alpha*order a cosine shaped
    % fade-out is used. Note, the fade-out of the original Tukey window is
    % defined such, that the windows is zero for m=order, if alpha~=0. The
    % window is modified such that the slope would yield zero for m=order+1
    alpha = conf.modal_window_parameter;
    m = ceil((1-alpha)*order):order;
    
    win = ones(1,order+1);
    win(m+1) = 0.5*(1 + cos(pi*(m-(1-alpha).*order)./(alpha*order+1)));
otherwise
    error('%s: unknown weighting type (%s)!',upper(mfilename),wtype);
end

% Inverse Circular Harmonics Transform
if nargout>1
    [varargout{1:nargout-1}] = inverse_cht([win(end:-1:2),win],ndtft);
end
