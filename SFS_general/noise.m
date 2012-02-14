function outsig = noise(samples,nsigs,noisetype)
%NOISE creates a white noise signal
%
%   Usage: outsig = noise(samples,nsigs,type)
%          outsig = noise(samples,nsigs)
%          outsig = noise(samples)
%
%   Input parameters:
%       samples     - length of the noise signal
%       nsigs       - number of noise signals, default: 1
%       type        - type of noise:
%                         'white' (default)
%                         'pink'
%                         'red'
%                         'brown'
%
%   Output parameters:
%       outsig  - samples x nsigs white noise signal, default: nsigs = 1. 
%
%   NOISE(samples,nsigs,type) generate a noise signal of type with a length of
%   samples and nsigs columns. The default value is nsigs = 1; the default value
%   for type is white.
%
%   NOTE: be sure you set the state of the random number generator in your
%   Matlab session, otherwise you will get the same noise every time you
%   restart Matlab.
%   The number generator is been set by:
%
%       randn('state',sum(100*clock));  % in Matlab < 7.7
%
%       s = RandStream.create('mt19937ar','seed',sum(100*clock)); 
%       RandStream.setDefaultStream(s);     % in Matlab >= 7.7
%
%   see also: click, irn

% AUTHOR: Hagen Wierstorf, Peter Soendergaard
% $LastChangedDate$
% $LastChangedRevision$
% $LastChangedBy$


%% ===== Checking input parameters =======================================
nargmin = 1;
nargmax = 3;
error(nargchk(nargmin,nargmax,nargin));
if nargin<3
    noisetype = 'white';
elseif nargin<2
    nsigs = 1;
end
isargpositivescalar(samples,nsigs)
isargchar(noisetype)


%% ===== Computation =====================================================
%
switch noisetype
    case 'white'
        % Generate white noise
        outsig = randn(samples,nsigs);

        % Generate white noise using FFT
        % Amplitude = 1
        %amplitude = ones(samples,1);
        % Generate random phase
        %phase = 2*pi*rand(samples,1);
        % Generate noise signal using fft
        %outsig = easyifft(amplitude,phase);

    case 'pink'
        % Handle trivial condition
        if samples==1
            outsig = ones(1,nsigs);
            return;
        end
        fmax = floor(samples/2)-1;
        f = (2:(fmax+1)).';
        % 1/f amplitude factor
        a = 1./sqrt(f);
        % Random phase
        p = randn(fmax,nsigs) + i*randn(fmax,nsigs);
        sig = repmat(a,1,nsigs).*p;
        outsig = ifftreal([ones(1,nsigs); sig; ...
            1/(fmax+2)*ones(1,nsigs)],samples);

    case 'red' || 'brown'
        outsig = cumsum(randn(samples,nsigs));

    otherwise
        error('%s: unknown noise type %s.',upper(mfilename),noisetype)
end

% Scale output noise signal
outsig = norm_signal(outsig);
