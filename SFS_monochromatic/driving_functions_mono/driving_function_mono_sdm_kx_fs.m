function D = driving_function_mono_sdm_kx_fs(kx,xs,f,conf)
%DRIVING_FUNCTION_MONO_SDM_KX_FS returns the driving signal D for a focused source in
%SDM in the kx domain
%
%   Usage: D = driving_function_mono_sdm_kx_fs(kx,xs,f,conf)
%
%   Input parameters:
%       kx          - kx dimension [nx1]
%       nk          - position of focused source / m [1x3]
%       f           - frequency of the monochromatic source / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_SDM_KX_FS(kx,xs,f,conf) returns SDM driving signals
%   for the given secondary sources, the focused source direction and the
%   frequency f. The driving signal is calculated in the kx domain.
%
%   References:
%       S. Spors and J. Ahrens (2010) - "Reproduction of Focused Sources by the
%       Spectral Division Method", ISCCSP
%
%   See also: driving_function_mono_wfs, driving_function_imp_wfs_ps

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
isargmatrix(kx,xs);
isargpositivescalar(f);
isargstruct(conf);


%% ===== Configuration ==================================================
xref = conf.xref;
c = conf.c;
dimension = conf.dimension;
driving_functions = conf.driving_functions;
x0 = conf.secondary_sources.center;
withev = conf.sdm.withev;


%% ===== Computation ====================================================
% Calculate the driving function in time-frequency domain

% Frequency
omega = 2*pi*f;
% Indexes for evanescent contributions and propagating part of the wave field
idxpr = (( abs(kx) <= (omega/c) ));
idxev = (( abs(kx) > (omega/c) ));
D = zeros(1,length(kx));

if strcmp('2D',dimension)

    % === 2-Dimensional ==================================================

    % Ensure 2D
    xs = xs(:,1:2);

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        to_be_implemented;
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2D focused source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

    % === 2.5-Dimensional ================================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        % D_25D(kx,w) = e^(i kx xs) ...
        %                                   ____________
        %                         H0^(2)( \|(w/c)^2-kx^2 |yref-ys| )
        %                     / - --------------_-_-_-_-_-_---------, |kx|<|w/c|
        %                     |      H0^(2)( \|(w/c)^2-kx^2 yref )
        %                    <        ____________
        %                     | K0( \|kx^2-(w/c)^2 |yref-ys| )
        %                     \ ----------_-_-_-_-_-_---------,       |kx|>|w/c|
        %                          K0( \|kx^2-(w/c)^2 yref )
        %
        % See Spors and Ahrens (2010), eq.(7)
        %
        D(idxpr) =  exp(1i*kx(idxpr)*xs(1)) .* ...
            besselh(0,2,sqrt( (omega/c)^2 - kx(idxpr).^2 )*abs(xref(2)-xs(2))) ./ ...
            besselh(0,2,sqrt( (omega/c)^2 - kx(idxpr).^2 )*abs(xref(2)-x0(2)));
        if(withev)
            D(idxev) =  exp(1i*kx(idxev)*xs(1)) .* ...
                besselk(0,sqrt(kx(idxev).^2 - (omega/c).^2)*abs(xref(2)-xs(2))) ./ ...
                besselk(0,sqrt(kx(idxev).^2 - (omega/c).^2)*abs(xref(2)-x0(2)));
        end

    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2.5D focused source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('3D',dimension)

    % === 3-Dimensional ==================================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        to_be_implemented;
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 3D focused source.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
