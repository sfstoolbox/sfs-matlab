function irs = new_irs()
%NEW_IRS creates an empty irs struct in the right format
%
%   Usage: irs = new_irs()
%
%   Output options
%       irs - irs struct in the desired format 
%
%   NEW_IRS() creates a irs struct in the format we have defined for the IR data
%   sets. This function creates a reference implementation of the format. So
%   don't change this function!
%
%   See also: check_irs, IR_format.txt

% AUTHOR: Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 0;
nargmax = 0;
error(nargchk(nargmin,nargmax,nargin));


%% ===== Computation =====================================================
% create our reference irs struct
irs = struct;
irs.description = 'Reference implementation of the irs struct.';
irs.head = 'dummy';             % Used dummy head
irs.ears = 'dummy';             % Used dummy head ears
irs.room = 'dummy';             % Used room
irs.source = 'dummy';           % Used loudsoeaker
irs.distance = 1;               % Distance between head and source. NOTE: this
                                %>has to be norm(head_position-source_position)
irs.fs = 44100;                 % Sampling rate
irs.head_position = [0 0 0];    % Position of head
irs.head_reference = [0 1 0];   % Position to which the head is pointing.
                                %>head_direction = head_reference-head_position
irs.source_position = [0 1 0];  % Position of loudspeaker source
irs.source_reference = [0 0 0]; % Position to which the source is pointing.
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
