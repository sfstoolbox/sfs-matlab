% wfs prefilter 3D (impulse response) result has to be applied to HRIR's over convolution theorem 
function hpre = wfs_prefilter3d(conf)


%% ===== Configuration ==================================================
fs = conf.fs;               % Sampling rate
flow = conf.wfs.hpreflow;   % Lower frequency limit of preequalization
                            % filter (= frequency when subwoofer is active)
fhigh = conf.wfs.hprefhigh; % Upper frequency limit of preequalization
                            % filter (= aliasing frequency of system)

%% ===== Variables ======================================================
Nfilt=128;
% Frequency axis
f = linspace(0,fs/2,fs/10);
% Find indices for frequencies in f smaller and nearest to fhigh and flow
idxfhigh = max(find(f<fhigh));
idxflow = max(find(f<flow));
% Initialize response
H = ones(1,length(f));


%% ===== Computation ====================================================

% Desired response for 3D
% Apply sqrt(2*pi*f)/sqrt(2*pi*fhigh) filter for idxf < idxfhigh
H(1:idxfhigh) = (2*pi*f(1:idxfhigh))./(2*pi*fhigh);
% % Set the response for idxf < idxflow to the value at idxflow
H(1:idxflow) = H(idxflow)*ones(1,idxflow);

% Compute filter
hpre = firls(Nfilt,2*f/fs,H);

% Truncate length to power of 2
hpre = hpre(1:end-1);


end
