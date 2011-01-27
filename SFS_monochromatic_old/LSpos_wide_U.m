% Loudspekeaker positions of U shaped flat panel array
% 
%
% ============================================================
% S.Spors / 19.11.03
% ============================================================
function LSpos = LSpos_wide_U(center_x,center_y)

LSdist = 0.171;          % exiter distance

LSpos = zeros(2,32);



LSpos(1,:) = -16*LSdist + center_x + linspace(0,32*LSdist,32);
LSpos(2,:) = center_y * ones(1,32);



