function test_delayline()
%TEST_DELAYLINE evaluates the delayline implementation
%
%   Usage: test_delayline()
%
%   TEST_DELAYLINE() tests the implementation of the delayline. For different
%   types (conf.usefracdelay) you have to manual edit the code below at the
%   moment.

%*****************************************************************************
% Copyright (c) 2010-2014 Quality & Usability Lab, together with             *
%                         Assessment of IP-based Applications                *
%                         Telekom Innovation Laboratories, TU Berlin         *
%                         Ernst-Reuter-Platz 7, 10587 Berlin, Germany        *
%                                                                            *
% Copyright (c) 2013-2014 Institut fuer Nachrichtentechnik                   *
%                         Universitaet Rostock                               *
%                         Richard-Wagner-Strasse 31, 18119 Rostock           *
%                                                                            *
% This file is part of the Sound Field Synthesis-Toolbox (SFS).              *
%                                                                            *
% The SFS is free software:  you can redistribute it and/or modify it  under *
% the terms of the  GNU  General  Public  License  as published by the  Free *
% Software Foundation, either version 3 of the License,  or (at your option) *
% any later version.                                                         *
%                                                                            *
% The SFS is distributed in the hope that it will be useful, but WITHOUT ANY *
% WARRANTY;  without even the implied warranty of MERCHANTABILITY or FITNESS *
% FOR A PARTICULAR PURPOSE.                                                  *
% See the GNU General Public License for more details.                       *
%                                                                            *
% You should  have received a copy  of the GNU General Public License  along *
% with this program.  If not, see <http://www.gnu.org/licenses/>.            *
%                                                                            *
% The SFS is a toolbox for Matlab/Octave to  simulate and  investigate sound *
% field  synthesis  methods  like  wave  field  synthesis  or  higher  order *
% ambisonics.                                                                *
%                                                                            *
% http://github.com/sfstoolbox/sfs                      sfstoolbox@gmail.com *
%*****************************************************************************


%% ===== Configuration ==================================================

% delays to evaluate
%dt=[-1 -0.9 -0.8 -0.7 -0.6 -0.5 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
%dt=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
dt=[-5 -2.5 0 2.5 5];
%dt=1*[0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1];
%dt=linspace(0,2,200);

% parameters for delayline
conf.usefracdelay = 1;
%conf.fracdelay_method = 'least_squares';
%conf.fracdelay_method = 'resample';
conf.fracdelay_method = 'interp1';

% length of input signal
L=256;

% create frequency axis
w = (0:1:(L-1))/L; 
wpi = w*pi;
wpi2=wpi(2:L);

% set up input signal
insig = zeros(L,1);
insig(L/2) = 1;


%% ===== Computation =====================================================
for n=1:length(dt)
    outsig(:,n) = delayline(insig,dt(n),1,conf);

    H(:,n) = freqz(outsig(:,n),1,wpi);
    magresp(:,n) = abs(H(:,n));
    uwphase(:,n)=-unwrap(angle(H(:,n)));
    
    phasdel(:,n) = uwphase(2:L,n)./wpi2';
end



%% ===== Plotting =====================================================
% setup legend and axis
t=1:L;
t=t-L/2;

% phase delay
figure;
plot(wpi2/pi,phasdel-(L/2)+1);
ylabel('phase delay');
xlabel('normalized frequency');
grid on;

% magnitude response
figure;
plot(wpi/pi,magresp);
ylabel('magnitude');
xlabel('normalized frequency');
grid on;

% impluse response
figure;
%plot(t(L/2-10:L/2+10),outsig(L/2-10:L/2+10,:));
imagesc(dt,t(L/2-50:L/2+50),db(abs(outsig(L/2-50:L/2+50,:))));
caxis([-100 10]);
ylabel('samples');
xlabel('delay');
turn_imagesc;
colorbar;
grid on;
