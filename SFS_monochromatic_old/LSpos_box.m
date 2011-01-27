% Loudspekeaker positions of box shaped array
% 
%
% ============================================================
% S.Spors / 09.12.03
% ============================================================
function [LSpos,LSdir] = LSpos_box(center_x,center_y,size,N_Ls)

LSdist = 4*size/N_Ls;          % exiter distance
N_Ls4 = N_Ls/4;


LSpos = zeros(2,N_Ls);
LSdir = zeros(1,N_Ls);


LSpos(1,1:N_Ls4) = center_x + linspace(0,size,N_Ls4) - N_Ls4/2*LSdist;
LSpos(2,1:N_Ls4) = center_y + ones(1,N_Ls4) * size/2 + LSdist;
% LSdir(1:N_Ls4) = pi;
LSdir(1:N_Ls4) = pi/2;


LSpos(1,N_Ls4+1:2*N_Ls4) = center_x + ones(1,N_Ls4) * size/2 + LSdist;
LSpos(2,N_Ls4+1:2*N_Ls4) = center_y - linspace(0,size,N_Ls4) + N_Ls4/2*LSdist;
LSdir(N_Ls4+1:2*N_Ls4) = pi;
% LSdir(N_Ls4+1:2*N_Ls4) = 3*pi/2;

LSpos(1,2*N_Ls4+1:3*N_Ls4) = center_x - linspace(0,size,N_Ls4) + N_Ls4/2*LSdist;
LSpos(2,2*N_Ls4+1:3*N_Ls4) = center_y - ones(1,N_Ls4) * size/2 - LSdist;
% LSdir(2*N_Ls4+1:3*N_Ls4) = 0;
LSdir(2*N_Ls4+1:3*N_Ls4) = 3*pi/2;


LSpos(1,3*N_Ls4+1:4*N_Ls4) = center_x - ones(1,N_Ls4) * size/2 - LSdist;
LSpos(2,3*N_Ls4+1:4*N_Ls4) = center_y + linspace(0,size,N_Ls4) - N_Ls4/2*LSdist;
LSdir(3*N_Ls4+1:4*N_Ls4) = 0;
% LSdir(3*N_Ls4+1:4*N_Ls4) = pi/2;


LSpos(2,:) = -LSpos(2,:);