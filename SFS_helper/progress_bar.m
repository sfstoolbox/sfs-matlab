function progress_bar(ii,nii,message)
%PROGRESS_BAR show the progress of an iteration
%
%   Usage: progress_bar(ii,nii,[message])
%
%   Input parameters:
%       ii      - current iteration
%       nii     - number of iterations
%       message - string printed before progress bar (default: 'Progress')

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


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 3;
narginchk(nargmin,nargmax);
isargpositivescalar(ii,nii);
if nargin<3
    message = 'Progress';
end


%% ===== Settings ========================================================
bar_width = 10;
bar_prefix = '[';
bar_suffix = ']';
bar_empty = ' ';
bar_full = '.';


%% ===== Generate the progress bar =======================================
% Format progress bar
progress = ii/nii;
nfull = floor(min(bar_width*progress,bar_width)); % number of full chars
nempty = bar_width - nfull;                       % number of empty chars
bar = repmat(bar_full,1,nfull);
empty = repmat(bar_empty,1,nempty);
% Print to screen
str = sprintf('%s: %s%s%s%s %3.0f%%%%', ...
              message, ...
              bar_prefix, ...
              bar, ...
              empty, ...
              bar_suffix, ...
              progress*100);
if ii==1
    fprintf(1,str);
else
    clear_line(length(str));
    fprintf(1,str);
end
if ii==nii
    fprintf(1,'\n');
end
% Octave didn't show the output directly in a function call, in order to do so
% it has explicitly flushed to stdout
if isoctave
    fflush(1);
end
end % end of function


%% ------ Subfunctions ---------------------------------------------------
function clear_line(str_length)
    % \r is not working under some Windows systems, therefore we use \b to
    % clear the line
    for ii=1:str_length
        fprintf(1,'\b');
    end
end
