function D = driving_function_mono_sdm_kx_ps(kx,xs,f,conf)
%DRIVING_FUNCTION_MONO_SDM_KX_PS returns the driving signal D for a point source in
%SDM in the kx domain
%
%   Usage: D = driving_function_mono_sdm_kx_ps(kx,xs,f,conf)
%
%   Input parameters:
%       kx          - kx dimension [nx1]
%       xs          - position of point source / m [1x3]
%       f           - frequency of the monochromatic source / Hz
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       D           - driving function signal [nx1]
%
%   DRIVING_FUNCTION_MONO_SDM_KX_PS(kx,xs,f,conf) returns SDM driving signals
%   for the given secondary sources, the virtual point source position and the
%   frequency f. The driving signal is calculated in the kx domain.
%
%   See also: driving_function_mono_sdm_kx

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
isargmatrix(kx);
isargxs(xs);
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
    xs = xs(1:2);

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        to_be_implemented;
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2D point source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

    % === 2.5-Dimensional ================================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        % D_25D(kx,w) = e^(i kx xs) ...
        %                                   ____________
        %                         H1^(2)( \|(w/c)^2-kx^2 |yref-ys| )
        %                     / - --------------_-_-_-_-_-_---------, |kx|<|w/c|
        %                     |      H1^(2)( \|(w/c)^2-kx^2 yref )
        %                    <        ____________
        %                     | K1( \|kx^2-(w/c)^2 |yref-ys| )
        %                     \ ----------_-_-_-_-_-_---------,       |kx|>|w/c|
        %                          K1( \|kx^2-(w/c)^2 yref )
        %
        % See http://sfstoolbox.org/#equation-D.sdm.ps.2.5D
        %
        D(idxpr) =  exp(1i*kx(idxpr)*xs(1)) .* ...
            besselh(1,2,sqrt( (omega/c)^2 - kx(idxpr).^2 )*abs(xref(2)-xs(2))) ./ ...
            besselh(1,2,sqrt( (omega/c)^2 - kx(idxpr).^2 )*abs(xref(2)-x0(2)));
        if(withev)
            D(idxev) =  exp(1i*kx(idxev)*xs(1)) .* ...
                besselk(1,sqrt(kx(idxev).^2 - (omega/c).^2)*abs(xref(2)-xs(2))) ./ ...
                besselk(1,sqrt(kx(idxev).^2 - (omega/c).^2)*abs(xref(2)-x0(2)));
        end
        %
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 2.5D point source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('3D',dimension)

    % === 3-Dimensional ==================================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        to_be_implemented;
    otherwise
        error(['%s: %s, this type of driving function is not implemented ', ...
            'for a 3D point source.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
