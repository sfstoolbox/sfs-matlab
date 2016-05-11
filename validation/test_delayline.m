function test_delayline()
%TEST_DELAYLINE evaluates the delayline implementation
%
%   Usage: test_delayline()
%
%   TEST_DELAYLINE() tests the implementation of the delayline. For different
%   types (conf.usefracdelay) you have to manual edit the code below at the
%   moment.

%*****************************************************************************
% The MIT License (MIT)                                                      *
%                                                                            *
% Copyright (c) 2010-2016 SFS Toolbox Team                                   *
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


%% ===== Configuration ==================================================
% Delays to evaluate
dt=[-5 -2.5 0 2.5 5];
%dt=[-1 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
%dt=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
%dt=1*[0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1];
%dt=linspace(0,2,200);
% Parameters for delayline
conf.usefracdelay = true;
fracdelay_methods = { ...
    'least_squares'; ...
    'resample'; ...
    'interp1'; ...
};

conf.fracdelay_method = 'interp1';

% Length of input signal
L=256;
% Create frequency axis
w = (0:1:(L-1))/L; 
wpi = w*pi;
wpi2=wpi(2:L);
% Set up input signal
insig = zeros(L,1);
insig(L/2) = 1;


%% ===== Computation and Plotting ========================================
for method = fracdelay_methods'
    conf.fracdelay_method = method{:};

    % --- Computation ---
    % Test all given delays
    for n=1:length(dt)
        outsig(:,n) = delayline(insig,dt(n),1,conf);
        H(:,n) = freqz(outsig(:,n),1,wpi);
        magresp(:,n) = abs(H(:,n));
        uwphase(:,n)=-unwrap(angle(H(:,n)));
        phasdel(:,n) = uwphase(2:L,n)./wpi2';
    end

    % --- Plotting ---
    % setup legend and axis
    t=1:L;
    t=t-L/2;
    % Phase delay
    figure;
    plot(wpi2/pi,phasdel-(L/2)+1);
    title([method{:},' - phase delay']);
    ylabel('phase delay');
    xlabel('normalized frequency');
    grid on;
    % Magnitude response
    figure;
    plot(wpi/pi,magresp);
    title([method{:},' - magnitude response']);
    ylabel('magnitude');
    xlabel('normalized frequency');
    grid on;
    % Impluse response
    figure;
    %plot(t(L/2-10:L/2+10),outsig(L/2-10:L/2+10,:));
    imagesc(dt,t(L/2-50:L/2+50),db(abs(outsig(L/2-50:L/2+50,:))));
    title([method{:},' - impulse response']);
    caxis([-100 10]);
    ylabel('samples');
    xlabel('delay');
    set(gca,'XTick',dt)
    turn_imagesc;
    colorbar;
    grid on;

end
