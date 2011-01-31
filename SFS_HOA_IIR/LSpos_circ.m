% Loudspekeaker positions and normal vectors for a circular array
% 
%
% ============================================================
% S.Spors / 26.7.07
% ============================================================
function [LSpos,LSdir] = LSpos_circ(center_x,center_y,size,N_LS)




LSpos = zeros(2,N_LS);
LSdir = zeros(1,N_LS);
phi = linspace(0,(1-1/N_LS)*2*pi,N_LS);

LSpos(1,:) = center_x+size*cos(phi);
LSpos(2,:) = center_y+size*sin(phi);

LSdir = unwrap(phi-pi);