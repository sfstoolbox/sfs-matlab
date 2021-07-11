% check_unified_25D_wfs() script
% add new functionality to the toolbox:
% we add and check the unified 2.5D WFS framework from
%   [Fir17] Gergely Firtha, Peter Fiala, Frank Schultz, Sascha Spors (2017):
%   "Improved Referencing Schemes for 2.5D Wave Field Synthesis Driving
%   Functions." In: IEEE/ACM Trans Audio Speech Lang Process,
%   25(5):1117-1127, DOI 10.1109/TASLP.2017.2689245.
%
%   [Sch17] Frank Schultz, Gergely Firtha, Peter Fiala, Sascha Spors (2017):
%   "Wave Field Synthesis Driving Functions for Large-Scale Sound
%   Reinforcement Using Line Source Arrays." In: Proc. of 142nd Audio Eng.
%   Soc. Conv., Berlin, #9722.
%
% note that [Fir17] radiates into +y for a linear SSD
% the simulations here consider radiation into -y for a linear SSD
%
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

if 0
    clear all;
    close all;
    clc;
    fh1 = figure;
    fh2 = figure;
else
    clc;
end

%% linear SSD:
%###############
conf = SFS_config;
conf.usetapwin = 0;
conf.secondary_sources.size = 60;
conf.secondary_sources.geometry = 'line';
conf.secondary_sources.number = 2^11;
f = 1500;
w_c = 2*pi*f/conf.c;

if 0
%% plane wave with linear SSD, [Fig. 3a/4a, Fir17]
nPW = [cos(300*pi/180) sin(300*pi/180) 0];
nPW = nPW/norm(nPW);

x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,nPW,'pw');
%either
dx0 = 2*ones(size(x0,1),1);  %[Fig. 3a, Fir17]
%or
%dx0 = abs(2/sin(300*pi/180))*ones(size(x0,1),1);  %[Fig. 4a, Fir17]

%% other referencing schemes:
%[(25), Fir17] for xref = const
%this handling is equivalent with Spors Revisited WFS (26),(27)
xRef = [0,-2,0];
%dx0 = zeros(size(x0,1),1);
for n=1:size(x0,1)
%   dx0(n) =  vector_norm(xRef-x0(n,1:3),2); 
end
%%

[P,x,y,z,x0,xPCS] = sound_field_mono_unified_25d_wfs([-6 6],[-3 0],0,nPW,'pw',f,conf,dx0);
P_ref = exp(-1i*w_c*(nPW(1)*x + nPW(2)*y));
figure(fh1)
surf(x,y,20*log10(abs(P-P_ref))), hold on
for n=1:size(x0,1)
   plot3(x0(n,1),x0(n,2),x0(n,2)*0+1e6,'ok','MarkerFaceColor','k','MarkerSize',.5) 
   plot3(xPCS(n,1),xPCS(n,2),xPCS(n,2)*0+1e3,'ow','MarkerFaceColor','w','MarkerSize',.5) 
end
hold off
shading flat
axis equal
axis([-6 6 -3 0])
view([0 90])
dBMax = 20;
dBMin = -35;
dBStep = 1;
N_cm = abs((dBMax-dBMin)/dBStep);
cm = (moreland(N_cm));
colormap(cm);
cb = colorbar('east');
cb.Ticks = [dBMin:5:dBMax];
set(gca,'CLim',[dBMin dBMax])
ylabel(cb,'dB')
xlabel('x / m')
ylabel('y / m')
figure(fh2)
surf(x,y,real(P))
shading flat
axis equal
axis([-6 6 -3 0])
view([0 90])
dBMax = 1;
dBMin = -1;
dBStep = 1/100;
N_cm = abs((dBMax-dBMin)/dBStep);
cm = (moreland(N_cm));
colormap(cm);
cb = colorbar('east');
cb.Ticks = [dBMin:0.25:dBMax];
set(gca,'CLim',[dBMin dBMax])
ylabel(cb,'dB')
xlabel('x / m')
ylabel('y / m')
end

if 0
%% line source with linear SSD, [Fig. 3b/4b, Fir17]
xPS = [0 1 0];

x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,xPS,'ls');

%either
dx0 = 2*ones(size(x0,1),1);  %[Fig. 3b, Fir17]
%or
%dx0 = ones(size(x0,1),1);  %prep for [Fig. 4b, Fir17]
%for n=1:size(dx0,1)
%   dx0(n,1) = dx0(n,1)*2*vector_norm(x0(n,1:3)-xPS,2)/xPS(2);  %[Fig. 4b, Fir17]
%end

%% other referencing schemes:
%[(25), Fir17] for xref = const
%this handling is equivalent with Spors Revisited WFS (26)
xRef = [0,-2,0];
%dx0 = zeros(size(x0,1),1);
for n=1:size(x0,1)
%   dx0(n) =  vector_norm(xRef-x0(n,1:3),2); 
end
%%

[P,x,y,z,x0,xPCS] = sound_field_mono_unified_25d_wfs([-6 6],[-3 0],0,xPS,'ls',f,conf,dx0);
r = sqrt((x-xPS(1)).^2+(y-xPS(2)).^2);
P_ref = -1i/4*besselh(0,2,w_c*r);
figure(fh1)
surf(x,y,20*log10(abs(P-P_ref))), hold on
for n=1:size(x0,1)
   plot3(x0(n,1),x0(n,2),x0(n,2)*0+1e6,'ok','MarkerFaceColor','k','MarkerSize',.5) 
   plot3(xPCS(n,1),xPCS(n,2),xPCS(n,2)*0+1e3,'ow','MarkerFaceColor','w','MarkerSize',.5) 
end
hold off
shading flat
axis equal
axis([-6 6 -3 0])
view([0 90])
dBMax = 0;
dBMin = -70;
dBStep = 1;
N_cm = abs((dBMax-dBMin)/dBStep);
cm = (moreland(N_cm));
colormap(cm);
cb = colorbar('east');
cb.Ticks = [dBMin:10:dBMax];
set(gca,'CLim',[dBMin dBMax])
ylabel(cb,'dB')
xlabel('x / m')
ylabel('y / m')
figure(fh2)
surf(x,y,real(P))
shading flat
axis equal
axis([-6 6 -3 0])
view([0 90])
dBMax = 0.1;
dBMin = -0.1;
dBStep = 1/1000;
N_cm = abs((dBMax-dBMin)/dBStep);
cm = (moreland(N_cm));
colormap(cm);
cb = colorbar('east');
cb.Ticks = [dBMin:1/50:dBMax];
set(gca,'CLim',[dBMin dBMax])
ylabel(cb,'dB')
xlabel('x / m')
ylabel('y / m')
end



if 1
%% point source with linear SSD, [Fig. 7, Fir17]
xPS = [0 3 0];
    
x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,xPS,'ps');
    
%either dx0=const
dx0 = 1.5*ones(size(x0,1),1);  %[Fig. 7a, Fir17], this would be the
%Ahrens WFS handling where generally dx0 = const is used (Sec. 3.9.3 in his book)

%or parallel reference line (equivalent with Start Diss (3.16)/(3.17) and HF/Far approx. SDM solution)
%dx0 = zeros(size(x0,1),1);
yref = -1.5; %[Fig. 7b, Fir17], [eq. (42), Fir17]
for n=1:size(x0,1)
%   dx0(n) =  vector_norm(x0(n,1:3)-xPS,2) * (-yref/(-yref-(-xPS(2))));
end

%or %[Fig. 7c, Fir17], [eq. (45), Fir17]
Rref = xPS(2)+1.5;
%dx0 = zeros(size(x0,1),1);
for n=1:size(x0,1)
   r0 = vector_norm(x0(n,1:3)-xPS,2);
   %dx0(n) =   r0*(Rref-r0)/Rref;
end

%% other referencing schemes:
xRef = [0,-1.5,0];
%dx0 = zeros(size(x0,1),1);
for n=1:size(x0,1) %from SPA 1, Schultz Diss (2.137), Start (Diss 3.10/3.11) with xref = const
   %dx0(n) =  vector_norm(xRef-x0(n,1:3),2)*vector_norm(x0(n,1:3)-xPS,2) / ...
   %    (vector_norm(xRef-x0(n,1:3),2) + vector_norm(x0(n,1:3)-xPS,2));
end

%Sascha Revisited WFS->uses the referencing function dx0 of 2D sound fields
%which for the 3D field of the virtual point source yields not the correct
%amplitude at desired reference positions (the SPA I ignores the influence of z-dimension)
%dx0 = zeros(size(x0,1),1);
for n=1:size(x0,1) %Sascha WFS revisited
   %dx0(n) =  vector_norm(xRef-x0(n,1:3),2)*2/3;  %the 2/3 factor compensates the amplitude mismatch for the chosen example
end

if 0
%one might use:
for n=1:size(x0,1)
   xRef = x0(n,1:3);
   xRef(2) = -3; %fixing yRef, but xref varies with x0
   dx0(n) =  vector_norm(xRef-x0(n,1:3),2)*vector_norm(x0(n,1:3)-xPS,2) / ...
       (vector_norm(xRef-x0(n,1:3),2) + vector_norm(x0(n,1:3)-xPS,2));
end
end
%%

[P,x,y,z,x0,xPCS] = sound_field_mono_unified_25d_wfs([-6 6],[-4 0],0,xPS,'ps',f,conf,dx0);
P_ref = exp(-1i*w_c*sqrt((x-xPS(1)).^2+(y-xPS(2)).^2)) ./(4*pi*sqrt((x-xPS(1)).^2+(y-xPS(2)).^2));
%%
figure(fh1)
surf(x,y,20*log10(abs(P-P_ref))), hold on
for n=1:size(x0,1)
   plot3(x0(n,1),x0(n,2),x0(n,2)*0+1e6,'ok','MarkerFaceColor','k','MarkerSize',.5) 
   plot3(xPCS(n,1),xPCS(n,2),xPCS(n,2)*0+1e3,'ow','MarkerFaceColor','w','MarkerSize',.5) 
end
hold off
shading flat
axis equal
axis([-6 6 -4 0])
view([0 90])
dBMax = 10;
dBMin = -80;
dBStep = 1;
N_cm = abs((dBMax-dBMin)/dBStep);
cm = (moreland(N_cm));
cm = parula(N_cm);
colormap(cm);
cb = colorbar('east');
cb.Ticks = [dBMin:10:dBMax];
set(gca,'CLim',[dBMin dBMax])
ylabel(cb,'dB')
xlabel('x / m')
ylabel('y / m')
figure(fh2)
surf(x,y,real(P))
shading flat
axis equal
axis([-6 6 -4 0])
view([0 90])
dBMax = 0.1;
dBMin = -0.1;
dBStep = 1/1000;
N_cm = abs((dBMax-dBMin)/dBStep);
cm = (moreland(N_cm));
colormap(cm);
cb = colorbar('east');
cb.Ticks = [dBMin:1/50:dBMax];
set(gca,'CLim',[dBMin dBMax])
ylabel(cb,'dB')
xlabel('x / m')
ylabel('y / m')

end



%% circular SSD:
%###############
conf = SFS_config;
conf.usetapwin = 0;
conf.secondary_sources.size = 3;
conf.secondary_sources.geometry = 'circular';
conf.secondary_sources.number = 2^9;
f = 2000;
w_c = 2*pi*f/conf.c;

if 0
%% plane wave with circular SSD, [Fig. 9a, Fir17] 
nPW = [1 0 0];
nPW = nPW/norm(nPW);

x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,nPW,'pw');
dx0 = 0.75*ones(size(x0,1),1);  %[Fig. 9/10, Fir17]

%% other referencing schemes:
%[(25), Fir17] for xref = const
%this handling is equivalent with Spors Revisited WFS (26),(27)
%if xRef is in origin we end up with dx0 = SSD radius
xRef = [-0.75,0,0];
%dx0 = zeros(size(x0,1),1);
for n=1:size(x0,1)
%   dx0(n) =  vector_norm(xRef-x0(n,1:3),2);
end
%%

[P,x,y,z,x0,xPCS] = sound_field_mono_unified_25d_wfs([-2 2],[-2 2],0,nPW,'pw',f,conf,dx0);
P_ref = exp(-1i*w_c*(nPW(1)*x + nPW(2)*y));
figure(fh1)
surf(x,y,20*log10(abs(P-P_ref))), hold on
for n=1:size(x0,1)
   plot3(x0(n,1),x0(n,2),x0(n,2)*0+1e6,'ok','MarkerFaceColor','k','MarkerSize',.5) 
   plot3(xPCS(n,1),xPCS(n,2),xPCS(n,2)*0+1e3,'ow','MarkerFaceColor','w','MarkerSize',.5) 
end
hold off
shading flat
axis square
axis([-2 2 -2 2])
view([0 90])
dBMax = 20;
dBMin = -35;
dBStep = 1;
N_cm = abs((dBMax-dBMin)/dBStep);
cm = (moreland(N_cm));
colormap(cm);
cb = colorbar('east');
cb.Ticks = [dBMin:5:dBMax];
set(gca,'CLim',[dBMin dBMax])
ylabel(cb,'dB')
xlabel('x / m')
ylabel('y / m')
end

if 0
%% line source with circular SSD, [Fig. 9b, Fir17]
xPS = [-2 0 0];

x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,xPS,'ls');
dx0 = 0.75*ones(size(x0,1),1);  %[Fig. 9/10, Fir17]

%% other referencing schemes:
%[(25), Fir17] for xref = const
%this handling is equivalent with Spors Revisited WFS (26)
%if xRef is in origin we end up with dx0 = SSD radius
xRef = [-0.75,0,0];
%dx0 = zeros(size(x0,1),1);
for n=1:size(x0,1)
%   dx0(n) =  vector_norm(xRef-x0(n,1:3),2);
end
%%

[P,x,y,z,x0,xPCS] = sound_field_mono_unified_25d_wfs([-2 2],[-2 2],0,xPS,'ls',f,conf,dx0);
r = sqrt((x-xPS(1)).^2+(y-xPS(2)).^2);
P_ref = -1i/4*besselh(0,2,w_c*r);
surf(x,y,20*log10(abs(P-P_ref))), hold on
for n=1:size(x0,1)
   plot3(x0(n,1),x0(n,2),x0(n,2)*0+1e6,'ok','MarkerFaceColor','k','MarkerSize',.5)
   plot3(xPCS(n,1),xPCS(n,2),xPCS(n,2)*0+1e3,'ow','MarkerFaceColor','w','MarkerSize',.5)    
end
hold off
shading flat
axis square
axis([-2 2 -2 2])
view([0 90])
dBMax = -10;
dBMin = -70;
dBStep = 1;
N_cm = abs((dBMax-dBMin)/dBStep);
cm = (moreland(N_cm));
colormap(cm);
cb = colorbar('east');
cb.Ticks = [dBMin:5:dBMax];
set(gca,'CLim',[dBMin dBMax])
ylabel(cb,'dB')
xlabel('x / m')
ylabel('y / m')
end

if 0
%% point source with circular SSD, [Fig. 10b, Fir17]
xPS = [-3 0 0]; %

x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,xPS,'ps');
dx0 = 0.75*ones(size(x0,1),1);  %[Fig. 9/10, Fir17], this would be the
%Ahrens WFS handling where generally dx0 = const is used (Sec. 3.9.3 in his book)

%% other referencing schemes:
% - Schultz Diss RefPoint Fig. 2.7 = SPA I
% - this is NOT equivalent to Sascha WFS revisited due to the 3D sound field
%   characteristics
xRef = [0,0,0];
%dx0 = zeros(size(x0,1),1);
for n=1:size(x0,1) %Schultz Diss reference point Fig. 2.7 
    % = SPA I = Start 1997 (3.10/3.11) = [(31), Fir17] for xref = const
   %dx0(n) =  vector_norm(xRef-x0(n,1:3),2)*vector_norm(x0(n,1:3)-xPS,2) / ...
   %    (vector_norm(xRef-x0(n,1:3),2) + vector_norm(x0(n,1:3)-xPS,2));
end

%Sascha Revisited WFS->uses the referencing function dx0 of 2D sound fields
%which for the 3D field of the virtual point source yields not the correct
%amplitude at desired reference positions (the SPA I ignores the influence of z-dimension)
for n=1:size(x0,1) %Sascha WFS revisited
   %dx0(n) =  vector_norm(xRef-x0(n,1:3),2); %*10^(-6.0206/20); %factor 1/2 corrects this primary source mismatch for the here chosen example  
end %and we end up precisely with the Ahrens approach if! xRef is in the origin, dx0 the corresponds to the SSD radius then
%%

[P,x,y,z,x0,xPCS] = sound_field_mono_unified_25d_wfs([-2 2],[-2 2],0,xPS,'ps',f,conf,dx0);
P_ref = exp(-1i*w_c*sqrt((x-xPS(1)).^2+(y-xPS(2)).^2)) ./(4*pi*sqrt((x-xPS(1)).^2+(y-xPS(2)).^2));

figure(fh1)
surf(x,y,20*log10(abs(P-P_ref))), hold on
for n=1:size(x0,1)
   plot3(x0(n,1),x0(n,2),x0(n,2)*0+1e6,'ok','MarkerFaceColor','k','MarkerSize',.5)
   plot3(xPCS(n,1),xPCS(n,2),xPCS(n,2)*0+1e3,'ow','MarkerFaceColor','w','MarkerSize',.5)    
end
hold off
shading flat
axis square
axis([-2 2 -2 2])
view([0 90])
dBMax = -10;
dBMin = -70;
dBStep = 1;
N_cm = abs((dBMax-dBMin)/dBStep);
cm = (moreland(N_cm));
colormap(cm);
cb = colorbar('east');
cb.Ticks = [dBMin:5:dBMax];
set(gca,'CLim',[dBMin dBMax])
ylabel(cb,'dB')
xlabel('x / m')
ylabel('y / m')
end



if 0
%% plane wave with circular SSD, [Fig. 12a, Fir17] 
nPW = [1 0 0];
nPW = nPW/norm(nPW);

x0 = secondary_source_positions(conf);
x0 = secondary_source_selection(x0,nPW,'pw');

%[eq. (51), Fir17] for [Fig. 12a, Fir17]
beta = 180 - 180/pi*cart2pol(x0(:,1),x0(:,2));
RSSD = beta*0 + conf.secondary_sources.size/2 ;
Rref = beta*0 +1.;
dx0 = RSSD.*cos(beta*pi/180) - sqrt(Rref.^2-RSSD.^2.*(sin(beta*pi/180)).^2); %[Fig. 12a, Fir17]
%complex dx0 here

%% other referencing schemes:
% - Schultz Diss RefPoint Fig. 2.7 = SPA I
% - this is equivalent to Sascha WFS revisited due to the 2D sound field
%   characteristics
% - this is in this special case (with xRef in the origin) also equivalent to
%   Ahrens WFS, where generally dx0 = const is used (Sec. 3.9.3 in his book) 
xRef = [0,0,0];
%dx0 = zeros(size(x0,1),1);
for n=1:size(x0,1)
   %dx0(n) =  vector_norm(xRef-x0(n,1:3),2); %Spors = Schultz = [(25), Fir17] for xref = const
   %dx0(n) =  1.5; % = Ahrens only if xRef is in origin
end
%%

[P,x,y,z,x0,xPCS] = sound_field_mono_unified_25d_wfs([-2 2],[-2 2],0,nPW,'pw',f,conf,dx0);
P_ref = exp(-1i*w_c*(nPW(1)*x + nPW(2)*y));
figure(fh1)
surf(x,y,20*log10(abs(P-P_ref))), hold on
for n=1:size(x0,1)
   plot3(x0(n,1),x0(n,2),x0(n,2)*0+1e6,'ok','MarkerFaceColor','k','MarkerSize',.5) 
   plot3(xPCS(n,1),xPCS(n,2),xPCS(n,2)*0+1e3,'ow','MarkerFaceColor','w','MarkerSize',.5) 
end
hold off
shading flat
axis square
axis([-2 2 -2 2])
view([0 90])
dBMax = 20;
dBMin = -30;
dBStep = 1;
N_cm = abs((dBMax-dBMin)/dBStep);
cm = (moreland(N_cm));
colormap(cm);
cb = colorbar('east');
cb.Ticks = [dBMin:5:dBMax];
set(gca,'CLim',[dBMin dBMax])
ylabel(cb,'dB')
xlabel('x / m')
ylabel('y / m')
end
