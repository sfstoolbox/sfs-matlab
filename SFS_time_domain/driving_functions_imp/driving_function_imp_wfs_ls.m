function [delay,weight] = driving_function_imp_wfs_ls(x0,nx0,xs,conf)
%DRIVING_FUNCTION_IMP_WFS_LS calculates the WFS weighting and delaying for a
%line source as source model
%
%   Usage: [delay,weight] = driving_function_imp_wfs_ls(x0,nx0,xs,conf)
%
%   Input parameters:
%       x0      - position  of secondary sources / m [nx3]
%       nx0     - direction of secondary sources [nx3]
%        xs     - position and orientation of virtual line source / m [nx3]
%                 or [nx6]
%       conf    - configuration struct (see SFS_config)
%
%   Output parameters:
%       delay   - delay of the driving function / s
%       weight  - weight (amplitude) of the driving function
%
%   DRIVING_FUNCTION_IMP_WFS_LS(x0,nx0,xs,conf) returns delays and weights for
%   the WFS driving function for a line source as source model. If no
%   explicit line source orientation xs(:,4:6) is given [0 0 1] is assumed.
%
%   References:
%       H. Wierstorf, J. Ahrens, F. Winter, F. Schultz, S. Spors (2015) -
%       "Theory of Sound Field Synthesis"
%
%   See also: sound_field_imp, sound_field_imp_wfs, driving_function_mono_wfs_ls

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
isargmatrix(x0,nx0,xs);
isargstruct(conf);


%% ===== Configuration ==================================================
% Speed of sound
c = conf.c;
xref = conf.xref;
fs = conf.fs;
dimension = conf.dimension;
driving_functions = conf.driving_functions;


%% ===== Computation =====================================================
% Get the delay and weighting factors

[xs,nxs] = get_position_and_orientation_ls(xs,conf);

if strcmp('2D',dimension)

    % === 2-Dimensional ==================================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        % d using a line source as source model
        %                     ___
        %                    | 1   (x0-xs) nx0
        % d(x0,t) = h(t) * _ |--- ------------- delta(t-|x0-xs|/c)
        %                   \|2pi |x0-xs|^(3/2)
        %
        % See http://sfstoolbox.org/#equation-d.wfs.ls
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Delay and amplitude weight
        delay = 1./c .* r;
        weight = 1./(2.*pi) .* vector_product(x0-xs,nx0,2) ./ r.^(3./2);
        %
    otherwise
        error(['%s: %s, this type of driving function is not implemented', ...
            'for a point source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

    % === 2.5-Dimensional ================================================

    % Reference point
    xref = repmat(xref,[size(x0,1) 1]);
    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        to_be_implemented;
    otherwise
        error(['%s: %s, this type of driving function is not implemented', ...
            'for a 2.5D point source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('3D',dimension)

    % === 3-Dimensional ==================================================

    switch driving_functions
    case 'default'
        % --- SFS Toolbox ------------------------------------------------
        % d using a line source as source model
        %                     ___
        %                    | 1   v nx0
        % d(x0,t) = h(t) * - |--- ------------- delta(t-|v|/c)
        %                   \|2pi |v|^(3/2)
        %
        % where v = x0-xs - <x0-xs,nxs > nxs,
        % and |nxs| = 1.
        %
        % See http://sfstoolbox.org/#equation-d.wfs.ls
        % and http://sfstoolbox.org/#equation-v.ls
        %
        % v = (I - nxs'nxs)(x0-xs)
        % r = |v|
        nxs = nxs(1,:);
        v = (x0 - xs)*(eye(3) - nxs'*nxs);
        r = vector_norm(v,2);
        delay = 1/c .* r;
        weight = 1/(2*pi) .* vector_product(v,nx0,2) ./ r.^(3/2);
    otherwise
        error(['%s: %s, this type of driving function is not implemented', ...
            'for a line source.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
