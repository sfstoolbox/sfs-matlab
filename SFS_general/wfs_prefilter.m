function hpre = wfs_prefilter(conf)
%WFS_PREFILTER generates a pre-equalization filter for WFS
%   Usage: hpre = wfs_prefilter(conf)
%          hpre = wfs_prefilter()

%
%   Input parameters:
%       conf    - optional struct containing configuration variables (see
%                 SFS_config for default values)
%
%   Output:
%       hpre    - preequalizaztion filter to convolute with BRIR 
%
%   WFS_PREFILTER(conf) calculates a sqrt(j k) pre-equalization filter for Wave 
%   Field Synthesis (from conf.hpreflow to conf.hprefhigh, see SFS_config).
%
%   see also: SFS_config, wfs_brs

% AUTHOR: Sascha Spors, Hagen Wierstorf


%% ===== Checking of input  parameters ==================================
nargmin = 0;
nargmax = 1;
error(nargchk(nargmin,nargmax,nargin));
if nargin<nargmax
    conf = SFS_config;
else
    isargstruct({conf},{'conf'});
end


%% ===== Configuration ==================================================

fs = conf.fs;               % Sampling rate

flow = conf.hpreflow;       % Lower frequency limit of preequalization 
                            % filter (= frequency when subwoofer is active)    
fhigh = conf.hprefhigh;     % Upper frequency limit of preequalization 
                            % filter (= aliasing frequency of system)

useplot = conf.useplot;     % Plot results?


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


%% ===== Plot resulting filter characteristics ==========================

if(useplot)
    Hfilt = fftshift(fft(hpre));
    Hfilt = Hfilt(Nfilt/2+1:end);
    f2 = linspace(0,fs/2,length(Hfilt));

    figure
    plot(f,20*log10(abs(H)));
    hold on
    plot(f2,20*log10(abs(Hfilt)),'ro-');
    hold off

    grid on
    xlabel('frequency -> [Hz]');
    ylabel('magnitude response -> [dB]');
    legend('desired response','filter response','Location','SouthEast');
    axis([0 2*fhigh -20 2]);


    figure
    freqz(hpre,1,[],fs)
end


%% ===== Save filter coefficients =======================================

if(0)
    fname = sprintf('wfs_prefilter_%d_%d_%d.wav',flow,fhigh,fs);
    disp(['Wrote pre-equalization filter into file: ' fname]);
    wavwrite(hpre,fs,fname);
end
