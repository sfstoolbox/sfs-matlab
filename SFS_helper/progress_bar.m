function progress_bar(ii,nii)
%PROGRESS_BAR show a progress of an iteration
%
%   Usage: progress_bar(ii,nii)
%
%   Input parameters:
%       ii  - current iteration
%       nii - number of iterations
%
%   PROGRESS_BAR(ii,nii) displays the progress of a loop.
%

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


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargpositivescalar(ii,nii);


%% ===== Generate the progress bar =======================================
% \r is not working under some Windows systems, therefore we use \b to clear the
% line
% calculate percentage
perc = ii/nii * 100;
% get bar
if perc<10
    bar = '[.         ]';
elseif perc<20
    bar = '[..        ]';
elseif perc<30
    bar = '[...       ]';
elseif perc<40
    bar = '[....      ]';
elseif perc<50
    bar = '[.....     ]';
elseif perc<60
    bar = '[......    ]';
elseif perc<70
    bar = '[.......   ]';
elseif perc<80
    bar = '[........  ]';
elseif perc<90
    bar = '[......... ]';
else
    bar = '[..........]';
end
str = sprintf('Progress: %s %3.0f%%%%',bar,perc);
if ii==1
    fprintf(1,str);
else
    clear_line(length(str)-1);
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
    for ii=1:str_length
        fprintf(1,'\b');
    end
end
