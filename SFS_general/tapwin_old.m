function win = tapwin(L,conf)
% NOTE: this is the old compatibility version in order to reproduce the stimuli
% for AES128, AES129
%TAPWIN generate a tapering window for a linear WFS array
%   Usage: win = tapwin(L,conf)
%          win = tapwin(L)
%
%   Input parameters:
%       L       - length of the loudspeaker array (m)
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output parameters:
%       win     - tapering window (1xnLS)
%
%   TAPWIN(L,conf) generates a tapering window for a linear WFS loudspeaker
%   array with a length of L.
%
%   see also: wfs_brs
%

% AUTHOR: Hagen Wierstorf, Sascha Spors


%% ===== Checking of input  parameters ==================================

if nargchk(1,2,nargin)
    error('Wrong number of args. Usage: win = tapwin(L,conf)');
end

if ~isnumeric(L) || ~isscalar(L) || L<0
    error('%s: L has to be a positive scalar!',upper(mfilename));
end

if nargin<2
    useconfig = true;
elseif ~isstruct(conf)
    error('%s: conf has to be a struct.',upper(mfilename));
else
    useconfig = false;
end


%% ===== Configuration ==================================================

% Load default configuration values
if(useconfig)
    conf = SFS_config;
end

% Apply tapering window?
usetapwin = conf.usetapwin;
tapwinlen = conf.tapwinlen;
% Loudspeaker distance
LSdist = conf.LSdist;

% Number of loudspeaker
nLS = ceil(L/LSdist);


%% ===== Calculation ====================================================
%
%    ------------------------------------------------------------ 
%   |                                                            |
% _|                                                              |_
%
if(usetapwin)
    
    % Length of window (ca. 30%, note: this will be splitted to the right 
    % and left site of the array)
    lenwin = ceil(tapwinlen*nLS)+1;
    % Create a hanning window with length lenwin
    %    -
    %   | |
    % _|   |_
    hannwin = hann(lenwin).^2;
    % Create tapering window
    win = [hannwin(1:ceil(end/2))' ...
           ones(1,nLS-lenwin) ...
           hannwin(ceil(end/2):end)'];
else
    % If you want to use no tapering window:
    win=ones(1,nLS);
end
