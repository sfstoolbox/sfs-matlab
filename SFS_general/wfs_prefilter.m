function hpre = wfs_prefilter(conf)
%WFS_PREFILTER creates a pre-equalization filter for WFS
%
%   Usage: hpre = wfs_prefilter(conf)
%          hpre = wfs_prefilter()
%
%
%   Input parameters:
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output:
%       hpre    - pre-equalization filter
%
%   WFS_PREFILTER(conf) calculates a sqrt(j k) pre-equalization filter for
%   Wave  Field Synthesis (from conf.hpreflow to conf.hprefhigh, see SFS_config).
%
%   see also: wfs_preequalization, SFS_config, brs_wfs_25d
%

% AUTHOR: Sascha Spors, Hagen Wierstorf
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$



%% ===== Checking of input  parameters ==================================
nargmin = 0;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct(conf);
end


%% ===== Configuration ==================================================
fs = conf.fs;               % Sampling rate
flow = conf.hpreflow;       % Lower frequency limit of preequalization 
                            % filter (= frequency when subwoofer is active)    
fhigh = conf.hprefhigh;     % Upper frequency limit of preequalization 
                            % filter (= aliasing frequency of system)


%% ===== Variables ======================================================

% Number of coefficients for filter
Nfilt=128;
% Frequency axis
f = linspace(0,fs/2,fs/10);
% Find indices for frequncies in f smaller and nearest to fhigh and flow
idxfhigh = max(find(f<fhigh));
idxflow = max(find(f<flow));
% Initialize response
H = ones(1,length(f));


%% ===== Computation ====================================================

% Desired response
% Apply sqrt(2*pi*f)/sqrt(2*pi*fhigh) filter for idxf < idxfhigh
H(1:idxfhigh) = sqrt(2*pi*f(1:idxfhigh))./sqrt(2*pi*fhigh);
% Set the response for idxf < idxflow to the value at idxflow
H(1:idxflow) = H(idxflow)*ones(1,idxflow);

% Compute filter
hpre = firls(Nfilt,2*f/fs,H);

% Truncate length to power of 2
hpre = hpre(1:end-1);
