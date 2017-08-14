function win = tapering_window(x0,conf)
%TAPERING_WINDOW generates a tapering window for a loudspeaker array
%
%   Usage: win = tapering_window(x0,conf)
%
%   Input parameters:
%       x0          - secondary sources / m
%       conf        - configuration struct (see SFS_config)
%
%   Output parameters:
%       win     - tapering window [nlsx1]
%
%   TAPERING_WINDOW(x0,conf) generates a tapering window for a secondary source
%   distribution given by x0. The window is created from a squared Hann window.
%   The strength of the tapering is controlled by the conf.tapwinlen setting.
%   If the secondary source distribution has some gaps, every joint part gets
%   its own tapering.
%
%   Tapering windows reduce diffraction in the synthesis sound field that
%   results from truncated secondary source distributions, see Sect. 3.2 of
%   Wierstorf (2014). The diffraction part of the sound field can be described
%   by cylindrical waves originating from the edges of the distribution, see
%   Sect. 8.3.2 in Born, Wolf (1999). Therefore a good method to reduce
%   truncation artifacts in the sound field is to fade out the amplitude of the
%   driving function at the edges of the array.
%
%   See also: secondary_source_position, sound_field_mono_wfs, hann
%
%
%   References:
%
%       Born, Wolf (1999) - "Principles of Optics", Cambridge University Press,
%       7th edition.
%
%       Wierstorf (2014) - "Perceptual Assessment of Sound Field Synthesis",
%       TU Berlin, https://doi.org/10.14279/depositonce-4310

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
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargsecondarysource(x0);
isargstruct(conf);


%% ===== Configuration ==================================================
usetapwin = conf.usetapwin;
tapwinlen = conf.tapwinlen;
geometry = conf.secondary_sources.geometry;


%% ===== Calculation =====================================================
% Number of secondary sources
nls = size(x0,1);
% Standard window with equal weights
win = ones(1,nls);
% FIXME: at the moment the tapering window is not working for spherical arrays,
% see https://github.com/sfstoolbox/sfs-matlab/issues/21
if usetapwin && nls>2 && ...
   ~(strcmp('sphere',geometry) || strcmp('spherical',geometry))
    % Get the mean distance between secondary sources and the smallest distance
    % to neighbour source for every secondary source. Due to long computing time
    % for really large secondary source numbers, the distance is approximated by
    % the first 100 (or less) sources. If you don't want this behavior, change
    % the following command to dx0 = secondary_source_distance(x0,0);
    dx0 = secondary_source_distance(x0,1);
    % Use only positions
    x0 = x0(:,1:3);
    % === Find the edges of the array ===
    edges = [];
    for ii=1:nls-1
        if norm(x0(ii,:)-x0(ii+1,:))>2*dx0
            edges = [edges; ii; ii+1];
        end
    end
    if norm(x0(end,:)-x0(1,:))>2*dx0
        edges = [edges; nls; 1];
    end
    % === Apply tapering window at all edges ===
    % If we have any edges in our array apply a tapering window for every array
    % part, consisting of two edges
    if ~isempty(edges)
        if edges(end)==1
            % First and last entry of secondary source is an edge
            edges = circshift(edges,[1,0]);
            start_idx = 1;
        else
            % First and last entry of secondary source is not an edge
            part_nls = edges(1) + nls-edges(end)+1;
            part_win = part_hann_win(part_nls,tapwinlen);
            win(1:edges(1)) = part_win(end-edges(1)+1:end);
            win(edges(end):end) = part_win(1:end-edges(1));
            start_idx = 2;
        end
        % Generate tapwin for every array part within the x0 vector
        for ii=start_idx:2:length(edges)-1
            part_nls = edges(ii+1)-edges(ii)+1;
            win(edges(ii):edges(ii+1)) = part_hann_win(part_nls,tapwinlen);
        end
    end
end

% Ensure the window will be a column vector
win = column_vector(win);

end

function [win] = part_hann_win(nls,tapwinlen)
    % Length of window (given by the value of tapwinlen). The window will be
    % splitted to both sides of the loudspeaker array.
    lenwin = round(tapwinlen*nls)+2;
    %
    % If we have less than four secondary sources, the usage of a tapering
    % window is no longer desired.
    if lenwin<4
        win = ones(1,nls);
    else
        % Create a squared Hann window with length lenwin
        win = hann_window(floor((lenwin-2)/2),floor((lenwin-2)/2),nls).^2;
    end
end
