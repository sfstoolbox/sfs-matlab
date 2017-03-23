function D = driving_function_mono_sdm_kx_pw(kx,nk,f,conf)
%DRIVING_FUNCTION_MONO_SDM_KX_PW returns the driving signal D for a plane wave in
%SDM in the kx domain
%
%   Usage: D = driving_function_mono_sdm_kx_pw(kx,nk,f,conf)
%
%   Input parameters:
%       kx          - kx dimension [nx1]
%       nk          - direction of plane wave / m [1x3]
%       f           - frequency of the monochromatic source / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_SDM_KX_PW(kx,nk,f,conf) returns SDM driving signals
%   for the given secondary sources, the virtual plane wave direction and the
%   frequency f. The driving signal is calculated in the kx domain.
%
%   See also: driving_function_mono_sdm_kx, sound_field_mono_sdm_kx

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
nargmin = 4;
nargmax = 4;
narginchk(nargmin,nargmax);
isargmatrix(kx,nk);
isargpositivescalar(f);
isargstruct(conf);


%% ===== Configuration ==================================================
xref = conf.xref;
c = conf.c;
dimension = conf.dimension;
driving_functions = conf.driving_functions;


%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

% Frequency
omega = 2*pi*f;
D = zeros(1,length(kx));

if strcmp('2D',dimension)

    % === 2-Dimensional ==================================================

    % Ensure 2D
    nk = nk(:,1:2);

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        to_be_implemented;
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2D plane wave.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

    % === 2.5-Dimensional ================================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        % D_2.5D using a plane wave as source model
        %
        %                   e^(-i w/c nky*xrefy)
        % D_2.5D(x0,w) = 4i ----------------------
        %                     (2) /w          \
        %                    H0  | - nky*xrefy |
        %                         \c          /
        %
        % See http://sfstoolbox.org/#equation-D.sdm.pw.2.5D
        %
        idx = find(kx>=omega/c*nk(:,1),1,'first');
        D(idx) = 4*1i*exp(-1i*omega/c*nk(2).*xref(2)) / ...
            besselh(0,2,omega/c.*nk(2).*xref(2));
        %
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2.5D plane wave.'],upper(mfilename),driving_functions);
    end


elseif strcmp('3D',dimension)

    % === 3-Dimensional ==================================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        to_be_implemented;
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 3D plane wave.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
