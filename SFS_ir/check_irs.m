function check_irs(irs)
%CHECK_IRS checks if the given irs is stored in the right format
%
%   Usage: check_irs(irs)
%
%   Input options
%       irs - irs struct
%
%   CHECK_IRS(irs) checks if the given irs is stored in our own format.
%   For format details have a look at the IR_format.txt file.
%
%   See also: new_irs, IR_format.txt

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
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);
isargstruct(irs);


%% ===== Format checking =================================================
% Check for the right field entries in the given irs struct
% Get a reference implementation of a irs and optional fields of this reference
% irs
[ref_irs,opt_fields] = new_irs();
% Get fields of reference implementation
ref_fields = fieldnames(ref_irs);
% Get the fields for the given irs
fields = fieldnames(irs);
idx = [];
for ii = 1:length(ref_fields)
    % Check if all needed fields are present
    if ~isfield(irs,ref_fields{ii}) && ~strcmp(opt_fields,ref_fields{ii})
        error('%s: The given irs misses the field: %s!',...
            upper(mfilename),ref_fields{ii});
    end
    % Get unneeded optional fields
    if ~isfield(irs,ref_fields{ii}) && strcmp(opt_fields,ref_fields{ii})
        idx = [idx ii];
    end
end
% Remove unneeded optional fields
ref_fields(idx) = [];

% Check if the order of the fields is standard conform to new_irs()
for ii = 1:length(ref_fields)
    if ~strcmp(fields{ii},ref_fields{ii})
        warning('SFS:irs_fields_order',...
                ['%s: the order of fields is not standard conform. ',...
                 'please use order_irs_fields(irs).'],upper(mfilename));
        break;
    end
end

% Check for right measurement angle types
if ~isnumeric(irs.head_azimuth) || ~isvector(irs.head_azimuth)
    error('%s: head_azimuth needs to be a vector.',upper(mfilename));
elseif ~isnumeric(irs.head_elevation) || ~isvector(irs.head_elevation)
    error('%s: head_elevation needs to be a vector.',upper(mfilename));
elseif ~isnumeric(irs.torso_azimuth) || ~isvector(irs.torso_azimuth)
    error('%s: torso_azimuth needs to be a vector.',upper(mfilename));
elseif ~isnumeric(irs.torso_elevation) || ~isvector(irs.torso_elevation)
    error('%s: torso_elevation needs to be a vector.',upper(mfilename));
end

% Check for the right number of entries for the signals and apparent angles
if size(irs.left,2)~=size(irs.right,2)
    error(['%s: the number of entries for the left ear signal is not ',...
           'consistent with the number of entries in the right channel.'],...
        upper(mfilename));
elseif size(irs.left,2)~=length(irs.apparent_azimuth)
    error(['%s: the number of entries for the left ear signal is not ',...
           'consistent with the number of entries for the apparent_azimuth.'],...
        upper(mfilename));
elseif size(irs.left,2)~=length(irs.apparent_elevation)
    error(['%s: the number of entries for the left ear signal is not ',...
           'consistent with the number of entries in the apparent_elevation.'],...
        upper(mfilename));
end

% Check for the right sizes of the entries for the positions
if all(~isnumeric(irs.head_position)) || size(irs.head_position,1)~=3
    error('%s: head_position needs to be a 3x1 vector.',upper(mfilename));
elseif all(~isnumeric(irs.head_reference)) || size(irs.head_reference,1)~=3 
     error('%s: head_reference needs to be a 3x1 vector.',upper(mfilename));
elseif all(~isnumeric(irs.source_position)) || size(irs.source_position,1)~=3
     error('%s: source_position needs to be a 3xn vector.',upper(mfilename));
elseif all(~isnumeric(irs.source_reference)) || size(irs.source_reference,1)~=3
    error('%s: source_reference needs to be a 3xn vector.',upper(mfilename));
%elseif isfield(irs,'room_corners')
%    if size(irs.room_corners,1)~=3
%        error('%s: room_corners needs to be a 3xn vector.',upper(mfilename));
%    end
end

% Check sampling rate
if ~isnumeric(irs.fs) || irs.fs<=0
    error('%s: fs needs to be a positive number.',upper(mfilename));
end

% Check distance
if isscalar(irs.distance)
    if ~isnumeric(irs.distance) || ...
        irs.distance-norm(irs.head_position-irs.source_position)>0.0001
        error('%s: distance has to be norm(head_position-source_position).',...
            upper(mfilename));
    end
else
    for ii = 1:length(irs.distance)
        if ~isnumeric(irs.distance(ii)) || ...
                irs.distance(ii) - ...
                norm(irs.head_position - irs.source_position(:,ii))>0.0001
            error(...
                '%s: distance has to be norm(head_position-source_position).',...
                upper(mfilename));
        end
    end
end

% Check string entries
if ~ischar(irs.description)
    error('%s: description needs to be a string.',upper(mfilename));
elseif ~ischar(irs.source)
    error('%s: loudspeaker needs to be a string.',upper(mfilename));
elseif ~ischar(irs.room)
    error('%s: room needs to be a string.',upper(mfilename));
elseif ~ischar(irs.head)
    error('%s: head needs to be a string.',upper(mfilename));
elseif ~ischar(irs.ears)
    error('%s: ears needs to be a string.',upper(mfilename));
end
