function [win,Win,Phi] = modal_weighting(order,ndtft,conf)
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
%       win         - the window w_n in the discrete domain (length = 2*order+1)
%       Win         - the inverse DTFT of w_n (length = ndtft)
%       Phi         - corresponding angle the inverse DTFT of w_n
%
%   MODAL_WEIGHTING(order,ndtft,conf) calculates a weighting window for the
%   modal band limitation applied in NFC-HOA. The window type is configured in
%   conf.nfchoa.wtype and its default setting is a simple rectangular window.
%   The other possibility is to use a Kaiser-Bessel window.
%
%   References:
%   	Kaiser, J., & Schafer, R. (1980) - "On the use of the I0-sinh window
%           for spectrum analysis", IEEE Transactions on Acoustics, Speech, and
%           Signal Processing
%       Van Trees, H. L. (2004) - "Optimum Array Processing", John Wiley & Sons.
%
%   See also: driving_function_imp_nfchoa, driving_function_mono_nfchoa

%*****************************************************************************
% Copyright (c) 2010-2016 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2016 Institut fuer Nachrichtentechnik                   *
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
wtype = conf.nfchoa.wtype;


%% ===== Computation =====================================================
switch wtype
    case 'rect'
        % === Rectangular Window =========================================
        win = ones(1,2*order+1);
    case {'kaiser', 'kaiser-bessel'}
        % === Kaiser-Bessel window =======================================
        % Approximation of the slepian window using modified bessel
        % function of zeroth order
        beta = conf.nfchoa.wparameter * pi;
        win = besseli(0,beta*sqrt(1-((-order:order)./order).^2)) ./ ...
              besseli(0,beta);
    otherwise
        error('%s: unknown weighting type (%s)!',upper(mfilename),wtype);
end

% Inverse DTFT
if nargout>1
    Win = ifft([win(order+1:end),zeros(1,order)],ndtft,'symmetric');
end
% Axis corresponding to inverse DTFT
if nargout>2
    Nphi = length(Win);
    Phi = 0:2*pi / Nphi:2*pi*(1-1/Nphi);
end
