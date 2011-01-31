clear all
clc
%==========================================================================
%% HOA (filter design) mainscript
% This script computates a (selectable) filter to get a reproduced sound
% field of a virtual acoustic scene with respect to the high order ambisonics approach.  
% The number of secondary sources (loudspeakers) is arbitrarily . The virtual acoustic scene
% can be choosen as plane wave (pw) or point source (ps).
% The filter can be choosen as 
%               (1) FIR filter 1 (computated with fir2 method)
%               (2) FIR filter 2 (computated with firintpolc method 
%                                 copyright 2009 by Hans W. Schuessler and Guenter F. Dehner)
%               (3) IIR filter   (computated regarding the diploma thesis of Pomberger)
%==========================================================================
%% simulation parameters
tic
config.nLS = 56;            % number of loudspeakers
config.al_pw = pi/2;        % incidence angle of plane wave
config.xps = [0 2];         % position of point/line source
config.xref = [0 0];        % reference position for 2.5D
config.c = 343;             % speed of sound
config.u = 300;             % number of sampling points

% frequency monofrequent
config.f = 1000;           
config.k = (2*pi*config.f)./config.c; 

% frequency broadband
config.f2 = linspace(1, 20000, config.u).';
config.k2 = ( 2.*pi ).*( config.f2 )./ (config.c);   

% radius of circular array / circular array
config.R = 1.50;  
[x0,n0] = LSpos_circ(0,0,config.R,config.nLS);
config.x0 = x0;
config.n0 = n0;

% filter parameters
config.N = 600;                          % filter order (for FIR filter)
config.N2 = 1024;                        % length of impulse (for IIR filter)
config.cLS = 14;                         % chosen LS for comparision
config.f3 = linspace(0, 1, config.u);    % discret sampling points regarding f2
config.n = 1024;                         % filter samples
config.fs = 44100;                       % sampling frequency
config.g = ones(length( config.f3 ), 1);% weighting factors
%==========================================================================
%% Comparison of Hankelfunction 2nd kind recursive and non recursive
% config.p = 1;   % choosen order for comparison 
% [w,config] = Bessel(config);
% [x,config] = hankel_rekursiv_real_imag(config);
%==========================================================================
%% calculation of HOA driving functions, broadband and monofrequent ps / pw
% with recursive calculation of Hankelfunctions
%==========================================================================

% [ps_b_r, ps_m_r] = HOA_driving_signal_broadband_mono_ps_rek(config);
% plot_frequency_response_ps_r;

% [pw_b_r, pw_m_r] = HOA_driving_signal_broadband_mono_pw_rek(config);
% plot_frequency_response_pw_r;

%==========================================================================
%% calculation of HOA driving functions, broadband and monofrequent ps / pw
% WITHOUT recursive calculation of Hankelfunctions
%==========================================================================

% [ps_b, ps_m] = HOA_driving_signal_broadband_mono_ps(config);
% plot_frequency_response_ps;

% [pw_b, pw_m] = HOA_driving_signal_broadband_mono_pw(config);
%  plot_frequency_response_pw;
 
% plot_comparison_of_frequency_responses_pw;
% plot_comparison_of_frequency_responses_ps;

%==========================================================================
%% filter design 
%==========================================================================
% chose ps or pw and rec or non rec
% E = ps_b;      % point source broadband , non recursive
% E = pw_b;      % plane wave broadband , non recursive
% E = ps_b_r;    % point source broadband recursive
% E = pw_b_r;    % plane wave broadband recursive
%==========================================================================
%% desired delay line for filter
%==========================================================================

% geometric
%delay = delay_line_fir_pw(config); %delay line for pw
%delay = delay_line_fir_ps(config); %delay line for ps

% group delay
% [delay] = group_delay(config, E);

%%=========================================================================
%
%% (1) fir2 - method (FIR filter)
%%==========================================================================

%%OUT returns desired coeeficients for every LS
% [OUT, X, F, Gd, F2] = HOA_filter_fir2(E, config,delay); 
% OUT = OUT ./ max(max(abs(OUT)));
% wavwrite(OUT,config.fs,'IRs\pw600fir2_groupDelay.wav');
%plot_frequency_response_fir2;

%==========================================================================

%% (2) firintpolc - method (config.n should be 597) (FIR filter)
%==========================================================================

%%OUT returns desired coeeficients for every LS
% [OUT, X, F, Gd, F2] = HOA_filter_firintpolc(E, config);
% OUT = OUT ./ max(max(abs(OUT)));
% wavwrite(OUT,config.fs,'IRs\ps597_firintpolc.wav');
% plot_frequency_response_firintpolc;

%==========================================================================
%% (3) IIR filter design
%==========================================================================
%[drivFuncCoeff] = IIR_laplace_pw(config);
%[drivFuncCoeff] = IIR_laplace_ps(config);
% plot_frequency_response_pw_neu;

%%OUT returns desired coeeficients for every LS
% [OUT] = IIR_filter_pw(config);
% OUT = OUT ./ max(max(abs(OUT)));
% wavwrite(OUT,config.fs,'IRs\iir_pw.wav');
% plot_freq_resp_and_imp_resp_of_iir_filter_pw;

%%OUT returns desired coeeficients for every LS
% [OUT] = IIR_filter_ps(config);
% OUT = OUT ./ max(max(abs(OUT)));
% wavwrite(OUT,config.fs,'IRs\iir_ps.wav');
% plot_freq_resp_and_imp_resp_of_iir_filter_ps;

%==========================================================================
%% plot reproduced wavefield from impulse response 
%% in frequency domain / time domain
%==========================================================================

% figure;
% plot_wave_field_from_imp_resp;
% figure;
% reproduced_field_time_domain_IR;

%==========================================================================
%% from here on monofrequent calculation of pw/ps necesarry
% you have to enable this in the section  "calculation of HOA driving
% functions, broadband and monofrequent ps / pw" then remove the return 
toc
return

%==========================================================================
%% compute reproduced wave field (monofrequent)
%==========================================================================

% simulation grid
x = linspace(-2,2,200);
y = linspace(-2,2,300);

P = zeros(length(x),length(y));
% choose ps or pw
% D = ps_m_r;
 D = pw_m_r;

idx=find(abs(D)~=0);
for n=idx
  P = P + D(n) .* point_source(x, y, config.x0(1,n), config.x0( 2,n ), config.k); % secondary point sources
 % P = P + D(n) .* line_source(x, y, config.x0(1,n), config.x0( 2,n ), config.k);  % secondary line sources
end
%==========================================================================
% plot reproduced wave field monofrequent
%==========================================================================
figure; 
plot_wavefield;
%==========================================================================
