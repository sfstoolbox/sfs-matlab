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
% Copyright (c) 2010-2013 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Deutsche Telekom Laboratories, TU Berlin           *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013      Institut fuer Nachrichtentechnik                   *
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
% http://dev.qu.tu-berlin.de/projects/sfs-toolbox       sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Checking of input parameters ====================================
nargmin = 2;
nargmax = 2;
narginchk(nargmin,nargmax);
isargpositivescalar(ii,nii);


%% ===== Generate the progress bar =======================================
% \r is not working under some Windows systems, therefore we use \b to clear the
% line
str = sprintf('Run %.0f/%.0f',ii,nii);
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
    for ii=1:str_length
        fprintf(1,'\b');
    end
end
