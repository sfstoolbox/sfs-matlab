function [irs,opt_fields] = new_irs()
%NEW_IRS creates an empty irs struct in the right format
%
%   Usage: irs = new_irs()
%
%   Output options
%       irs         - irs struct in the desired format
%       opt_fields  - cell array containing which of the fields in the returning
%                     irs struct are optional
%
%   NEW_IRS() creates a irs struct in the format we have defined for the IR data
%   sets. This function creates a reference implementation of the format. So
%   don't change this function!
%
%   See also: check_irs, IR_format.txt

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
nargmin = 0;
nargmax = 0;
narginchk(nargmin,nargmax);


%% ===== Computation =====================================================
% the following fields are only optional
opt_fields = {'room_corners'};
% create our reference irs struct
irs = struct;
irs.description = 'Reference implementation of the irs struct.';
irs.head = 'dummy';             % Used dummy head
irs.ears = 'dummy';             % Used dummy head ears
irs.room = 'dummy';             % Used room
irs.room_corners = [0 0 0]';    % corners of the used room
irs.source = 'dummy';           % Used loudsoeaker
irs.distance = 1;               % Distance between head and source. NOTE: this
                                %>has to be norm(head_position-source_position)
irs.fs = 44100;                 % Sampling rate
irs.head_position = [0 0 0]';   % Position of head
irs.head_reference = [0 1 0]';  % Position to which the head is pointing.
                                %>head_direction = head_reference-head_position
irs.source_position = [0 1 0]'; % Position of loudspeaker source
irs.source_reference = [0 0 0]';% Position to which the source is pointing.
                                %>source_direction =
                                %>source_reference-source_position
irs.head_azimuth = NaN;         % Head azimuth (NaN if no rotation took place)
irs.head_elevation = NaN;       % Head elevation (NaN if no rotation took place)
irs.torso_azimuth = NaN;        % Torso azimuth (NaN if no rotation took place)
irs.torso_elevation = NaN;      % Torso elevation (NaN if no rotation took
                                %>place)
irs.apparent_azimuth = [];      % Apparent azimuth of source in relation to
                                %>head direction
irs.apparent_elevation = [];    % Apparent elevation of source in relation to
                                %>head direction
irs.left = [];                  % Left ear signal
irs.right = [];                 % Right ear signal
