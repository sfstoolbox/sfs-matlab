function status = test_wfs_25d(modus)
%TEST_WFS_25D tests behavior of 2.5D WFS
%
%   Usage: status = test_wfs_25d(modus)
%
%   Input parameters:
%       modus   - 0: numerical
%                 1: visual
%
%   Output parameters:
%       status  - true or false

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


status = false;


%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);


%% ===== Configuration ===================================================
% Parameters
conf = SFS_config;
conf.secondary_sources.geometry = 'linear';
conf.secondary_sources.size = 20;
conf.secondary_sources.number = 512;
conf.secondary_sources.center = [0,6,0];
conf.xref = [0,0,0];
conf.usetapwin = true;
conf.tapwinlen = 0.2;

%
f = 2000;  % temporal frequency
positions = { [0,9,0], [0,-1,0], [0,3,0,0,-1,0] };  % source positions
sources = {'ps', 'pw', 'fs'};
gtsources = {'ps', 'pw', 'ps'};
X = [-2,2];
Y = [-2,2];
Z = 0;


%% ===== Main ============================================================
for idx=1:length(positions)
    xs = positions{idx};
    src = sources{idx};
    gt = gtsources{idx};

    if modus
        figure;
        ddx = 0;
    end

    for driving_functions = {'reference_point', 'reference_line'}

        conf.driving_functions = driving_functions{:};
        Pgt = sound_field_mono(X,Y,Z,[xs(1:3),0,-1,0,1],gt,1,f,conf);
        Pwfs = sound_field_mono_wfs(X,Y,Z,xs,src,f,conf);

        if modus
            subplot(2,2,2*ddx+1);
            imagesc(Y,X,real(Pwfs));
            title(sprintf('%s %s',src,driving_functions{:}),'Interpreter','none');
            set(gca,'YDir','normal');
            colorbar;

            subplot(2,2,2*ddx+2);
            imagesc(Y,X,real(db(1 - Pwfs./Pgt)));
            title(sprintf('%s %s',src,driving_functions{:}),'Interpreter','none');
            set(gca,'YDir','normal');
            colorbar;
            hold on;
            if strcmp('reference_point',conf.driving_functions)
                plot(conf.xref(1),conf.xref(2),'gx');
            else
                plot(conf.xref(1)+X,conf.xref([2,2]),'g--');
            end
            hold off;

            ddx= ddx+1;
        end
    end
end


status = true;
