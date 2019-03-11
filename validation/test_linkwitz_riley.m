function status = test_linkwitz_riley(modus)
%TEST_LINKWITZ_RILEY tests the design of the Linkwitz-Riley Filters
%
%   Usage: status = test_linkwitz_riley(modus)
%
%   Input parameters:
%       modus   - 0: numerical
%                 1: visual
%
%   Output parameters:
%       status  - true or false

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2019 SFS Toolbox Developers                             *
%                                                                            *
% Permission is hereby granted,  free of charge,  to any person  obtaining a *
% copy of this software and associated documentation files (the "Software"), *
% to deal in the Software without  restriction, including without limitation *
% the rights  to use, copy, modify, merge,  publish, distribute, sublicense, *
% and/or  sell copies of  the Software,  and to permit  persons to whom  the *
% Software is furnished to do so, subject to the following conditions:       *
%                                                                            *
% The above copyright notice and this permission notice shall be included in *
% all copies or substantial portions of the Software.                        *
%                                                                            *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
% IMPLIED, INCLUDING BUT  NOT LIMITED TO THE  WARRANTIES OF MERCHANTABILITY, *
% FITNESS  FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO EVENT SHALL *
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER *
% LIABILITY, WHETHER  IN AN  ACTION OF CONTRACT, TORT  OR OTHERWISE, ARISING *
% FROM,  OUT OF  OR IN  CONNECTION  WITH THE  SOFTWARE OR  THE USE  OR OTHER *
% DEALINGS IN THE SOFTWARE.                                                  *
%                                                                            *
% The SFS Toolbox  allows to simulate and  investigate sound field synthesis *
% methods like wave field synthesis or higher order ambisonics.              *
%                                                                            *
% https://sfs.readthedocs.io                            sfstoolbox@gmail.com *
%*****************************************************************************

status = false;

%% ===== Checking of input  parameters ===================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);

%% ===== Main ============================================================
fs = 44100;  % sampling frequency in Hz
fc = 1000;  % crossover frequency in Hz
f = logspace(-3,log10(fs/2),10000);  % frequeny-axis in Hz
omega = 2*pi*f;  % angular frequency

ftype = {'low' 'high' 'all'};  % filter types
Hs = zeros(3,length(f));  % s-domain frequency spectra
Hz = Hs;  % z-domain frequency spectra

for Nlr = [2 4 12]  % filter orders
  
  for tdx = 1:3
    % Laplace-Domain
    [zs,ps,ks] = linkwitz_riley(Nlr,2*pi*fc,ftype{tdx},'s');
    [soss,gs] = zp2sos(zs,ps,ks,'down','none');
    
    s = [-omega.^2; 1i*omega];
    s(3,:) = 1;
    
    if isoctave && strcmp(ftype{tdx}, 'low')
      % WORKAROUND for Octave bug (https://savannah.gnu.org/bugs/?51936)
      Hs(tdx,:) = gs.*prod( (soss(:,3:-1:1)*s)./(soss(:,4:6)*s),1);
    else
      Hs(tdx,:) = gs.*prod( (soss(:,1:3)*s)./(soss(:,4:6)*s),1);
    end
    
    % z-Domain
    [zz,pz,kz] = linkwitz_riley(Nlr,fc./fs*2,ftype{tdx},'z');
    [sosz,gz] = zp2sos(zz,pz,kz,'down','none');
    
    z = [exp(0.*omega./fs); exp(-1j.*omega./fs); exp(-2j.*omega./fs)];
    Hz(tdx,:) = gz.*prod( (sosz(:,1:3)*z)./(sosz(:,4:6)*z),1);
  end
  
  % plotting
  if modus
    figure;
    subplot(2,2,1);
    plot_amp(f, Hs);
    title(sprintf('Amplitude Spectrum, s-Domain, order=%d', Nlr));
    subplot(2,2,2);
    plot_phase(f, Hs);
    title(sprintf('Phase Spectrum, s-Domain, order=%d', Nlr));
    subplot(2,2,3);
    plot_amp(f, Hz);
    title(sprintf('Amplitude Spectrum, z-Domain, order=%d', Nlr));
    subplot(2,2,4);
    plot_phase(f, Hz);
    title(sprintf('Phase Spectrum, z-Domain, order=%d', Nlr));
  end
end

status = true;

end

function plot_amp(f, H)
semilogx( ...
  f, db(H(1,:)), 'b', ...
  f, db(H(2,:)), 'r', ...
  f, db(H(1,:)+H(2,:)), 'g.', ...
  f, db(H(3,:)), 'k--' ...
  );
legend('Lowpass', 'Highpass', 'Lowpass + Highpass', 'Allpass', 'Location', 'southwest');
xlabel('f / Hz');
ylabel('20 lg |H| / dB')
xlim([f(1), f(end)]);
ylim([-90, 10]);
end

function plot_phase(f, H)
semilogx( ...
  f, unwrap(angle(H(1,:))), 'b', ...
  f, unwrap(angle(H(2,:))), 'r', ...
  f, unwrap(angle(H(1,:)+H(2,:))), 'g.', ...
  f, unwrap(angle(H(3,:))), 'k--' ...
  );
legend('Lowpass', 'Highpass', 'Lowpass + Highpass', 'Allpass', 'Location', 'southwest');
xlabel('f / Hz');
ylabel('angle(H) / rad')
xlim([f(1), f(end)]);
end
