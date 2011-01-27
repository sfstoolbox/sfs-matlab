% Loudspekeaker positions for a linear array
% 
%
% ============================================================
% S.Spors / 09.12.03
% ============================================================
function [LSpos,LSdir] = LSpos_box(center_x,center_y,size,nLS)


x=linspace(-size/2,size/2,nLS);

LSpos(1,:) = center_x + linspace(-size/2,size/2,nLS);
LSpos(2,:) = center_y * ones(1,nLS);


LSdir = pi/2*ones(1,nLS);


