function status = test_delayline(modus)
%TEST_DELAYLINE evaluates the delayline implementation
%
%   Usage: status = test_delayline(modus)
%
%   Input parameters:
%       modus   - 0: numerical
%                 1: visual
%
%   Output parameters:
%       status  - true or false
%
%   TEST_DELAYLINE(modus) tests the implementation of the delayline. For
%   different types (conf.usefracdelay) you have to manual edit the code below
%   at the moment.

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2017 SFS Toolbox Developers                             *
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
% http://sfstoolbox.org                                 sfstoolbox@gmail.com *
%*****************************************************************************


status = false;


%% ===== Checking of input parameters ====================================
nargmin = 1;
nargmax = 1;
narginchk(nargmin,nargmax);


%% ===== Configuration ==================================================
% Delays to evaluate in samples (they will be transformed into seconds later)
dt = [-5 -2.5 0 2.5 5];
%dt=[-1 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
%dt=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
%dt=1*[0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1];
%dt=linspace(0,2,200);

% Length of input signal
L=256;
% Create frequency axis
w = (0:1:(L-1))/L;
wpi = w*pi;
wpi2=wpi(2:L);
% Set up input signal
insig = zeros(L,1);
insig(L/2) = 1;

conf.fs = 44100;
conf.delayline.filterorder = 4;  % / 1
conf.delayline.resamplingfactor = 4;  % / 1
conf.delayline.resamplingorder = 128; % / 1

%% ===== Computation and Plotting ========================================
for resampling = {'none', 'matlab', 'pm'}
    conf.delayline.resampling = resampling{:};
    for filter = {'lagrange', 'thiran', 'least_squares', 'integer'}
        conf.delayline.filter = filter{:};

        % --- Computation ---
        % Test all given delays
        for n=1:length(dt)
            [outsig(:,n), delay_offset] = delayline(insig,dt(n)/conf.fs,1,conf);
            H(:,n) = freqz(outsig(:,n),1,wpi);
            magresp(:,n) = abs(H(:,n));
            uwphase(:,n)=-unwrap(angle(H(:,n)));
            phasdel(:,n) = uwphase(2:L,n)./wpi2';
        end

        % --- Plotting ---
        if modus
            figure;
            % setup legend and axis
            t=1:L;
            t=t-L/2;
            % Phase delay
            subplot(2,2,1);
            plot(wpi2/pi,phasdel-(L/2)+1-delay_offset);
            title(['resample: ', resampling{:}, ', filter: ', filter{:},' - phase delay']);
            ylabel('phase delay');
            xlabel('normalized frequency');
            legend(num2str(dt.','%.1f'));
            grid on;
            % Magnitude response
            subplot(2,2,2);
            plot(wpi/pi,magresp);
            title(['resample: ', resampling{:}, ', filter: ', filter{:},' - magnitude response']);
            ylabel('magnitude');
            xlabel('normalized frequency');
            legend(num2str(dt.','%.1f'));
            grid on;
            % Impluse response
            subplot(2,2,3);
            imagesc(dt,t(L/2-50:L/2+50),db(abs(outsig(L/2-50:L/2+50,:))));
            title(['resample: ', resampling{:}, ', filter: ', filter{:},' - impulse response']);
            caxis([-100 10]);
            ylabel('samples');
            xlabel('delay / samples');
            set(gca,'XTick',dt)
            turn_imagesc;
            colorbar;
            grid on;
        end
    end
end


status = true;
