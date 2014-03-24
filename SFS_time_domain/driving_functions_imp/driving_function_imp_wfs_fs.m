function [delay,weight] = driving_function_imp_wfs_fs(x0,nx0,xs,conf)
%DRIVING_FUNCTION_IMP_WFS_FS calculates the WFS weighting and delaying for a
%focused source as source model
%
%   Usage: [delay,weight] = driving_function_imp_wfs_fs(x0,nx0,xs,[conf]);
%
%   Input parameters:
%       x0      - position  of secondary sources (m) [nx3]
%       nx0     - direction of secondary sources [nx3]
%       xs      - position of focused source [nx3]
%       conf    - optional configuration struct (see SFS_config)
%
%   Output parameters:
%       delay   - delay of the driving function (s)
%       weight  - weight (amplitude) of the driving function
%
%   DRIVING_FUNCTION_IMP_WFS_FS(x0,nx0,xs,conf) returns delays and weights for
%   the WFS driving function for a focused source as source model.
%
%   References:
%       H. Wierstorf (2014) - "Perceptual Assessment of Sound Field Synthesis",
%       PhD thesis, Tu Berlin
%
%   see also: sound_field_imp, sound_field_imp_wfs, driving_function_mono_wfs_fs

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
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


%% ===== Checking of input  parameters ==================================
nargmin = 3;
nargmax = 4;
narginchk(nargmin,nargmax);
isargmatrix(x0,nx0,xs);
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
% Speed of sound
c = conf.c;
xref = conf.xref;
fs = conf.fs;
dimension = conf.dimension;
driving_functions = conf.driving_functions;


%% ===== Computation =====================================================

% Get the delay and weighting factors
if strcmp('2D',dimension)

    % === 2-Dimensional ==================================================

    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        to_be_implemented;
        % d using a focused line source as source model
        %                     ___
        %                    | 1   (x0-xs) nx0
        % d(x0,t) = h(t) * - |--- ------------- delta(t+|x0-xs|/c)
        %                   \|2pi |x0-xs|^(3/2)
        %
        % this is a time reversered version of Wierstorf2014 p.26, (2.57)
        %
        % r = |x0-xs|
        r = vector_norm(x0-xs,2);
        % Delay and amplitude weight
        delay = -1/c .* r;
        weight = 1/(2*pi) .* vector_product(x0-xs,nx0,2) ./ r.^(3/2);
    else
        error(['%s: %s, this type of driving function is not implemented', ...
            'for a 2D focused source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('2.5D',dimension)

    % === 2.5-Dimensional ================================================

    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % Reference point
        xref = repmat(xref,[size(x0,1) 1]);
        % 2.5D correction factor
        %        ______________
        % g0 = \| 2pi |xref-x0|
        %
        g0 = sqrt(2*pi*vector_norm(xref-x0,2));
        %
        % d_2.5D using a line sink as source model
        %
        %                       -g0 (xs-x0) nx0
        % d_2.5D(x0,t) = h(t) * --- ------------- delta(t + |xs-x0|/c)
        %                       2pi |xs-x0|^(3/2)
        %                       
        % see Wierstorf (2014), p.27 (2.65)
        %
        % r = |xs-x0|
        r = vector_norm(xs-x0,2);
        % Delay and amplitude weight
        delay =  -1/c .* r;
        weight = -g0/(2*pi) .* vector_product(xs-x0,nx0,2) ./ r.^(3/2);
    else
        error(['%s: %s, this type of driving function is not implemented', ...
            'for a 2.5D focused source.'],upper(mfilename),driving_functions);
    end


elseif strcmp('3D',dimension)

    % === 3-Dimensional ==================================================

    if strcmp('default',driving_functions)
        % --- SFS Toolbox ------------------------------------------------
        % d_3D using a point sink as source model
        % 
        %                     -1  (xs-x0) nx0
        % d_3D(x0,t) = h(t) * --- ------------- delta(t + |xs-x0|/c)
        %                     2pi |xs-x0|^(3/2)
        %
        % see Wierstorf (2014), p.27 (2.64)
        %
        % r = |xs-x0|
        r = vector_norm(xs-x0,2);
        % Delay and amplitude weight
        delay = -1/c .* r;
        weight = -1/(2*pi) .* vector_product(xs-x0,nx0,2) ./ r.^(3/2);
    else
        error(['%s: %s, this type of driving function is not implemented', ...
            'for a 3D focused source.'],upper(mfilename),driving_functions);
    end

else
    error('%s: the dimension %s is unknown.',upper(mfilename),dimension);
end
