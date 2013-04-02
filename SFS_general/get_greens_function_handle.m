function h = get_greens_function_handle(src,domain)
%GET_GREENS_FUNCTION_HANDLE returns the handle to the desired Green's function
%
%   Usage: h = get_greens_function_handle(src,domain)
%
%   Input options:
%       src     - Green's function corresponds to this source type.
%                 Valid types are:
%                   'ps' - point source
%                   'ls' - line source
%                   'pw' - plane wave
%       domain  - Green's function for the time or frequency domain. This is
%                 indicated by setting domain to 't' or 'f'
%
%   Output options:
%       h       - function handle to the desired Green's function
%
%   GET_GREENS_FUNCTION_HANDLE(src,domain) returns a function handle to the
%   Green's function corresponding to the desired source type and domain.
%   For example:
%   h = get_greens_function_handle('ps','f') returns a handle to the function
%   point_source_mono(), which can then be called via h().
%
%   see also: wave_field_mono, wave_field_imp

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


%% ===== Computation =====================================================
% Greens functions for the time domain
if strcmp('t',domain) || strcmp('time',domain)
    if strcmp('ps',src)
        h = @point_source_imp;
    elseif strcmp('ls',src)
        h = @line_source_imp;
    elseif strcmp('pw',src)
        h = @plane_wave_imp;
    else
        error('%s: %s is not a valid source model',upper(mfilename),dimension);
    end
% Greens functions for the frequency domain    
elseif strcmp('f',domain) || strcmp('frequency',domain)
    if strcmp('ps',src)
        h = @point_source_mono;
    elseif strcmp('ls',src)
        h = @line_source_mono;
    elseif strcmp('pw',src)
        h = @plane_wave_mono;
    else
        error('%s: %s is not a valid source model',upper(mfilename),dimension);
    end
else
    error('%s: %s is not a valid domain, use "t" or "f".', ...
        upper(mfilename),domain);
end
