% omega = 2*pi*logspace(-1,4,1000);
% omega_vec = [-omega.^2; 1i*omega];
% omega_vec(3,:) = 1;
% 
% Nlr = 12;
% 
% [zap,pap,kap] = linkwitz_riley(Nlr,500,'all','s');  % LR Allpass filter
% [sosap,g] = zp2sos(zap,pap,kap,'down','none');  % SOS of Allpass
% Hap = g.*prod( (sosap(:,1:3)*omega_vec)./(sosap(:,4:6)*omega_vec),1);
% 
% [zhp,php,khp] = linkwitz_riley(Nlr,500,'high','s');  % LR Allpass filter
% [soshp,g] = zp2sos(zhp,php,khp,'down','none');  % SOS of Allpass
% Hhp = g.*prod( (soshp(:,1:3)*omega_vec)./(soshp(:,4:6)*omega_vec),1);
% 
% [zlp,plp,klp] = linkwitz_riley(Nlr,500,'low','s');  % LR Allpass filter
% [soslp,g] = zp2sos(zlp,plp,klp,'down','none');  % SOS of Allpass
% Hlp = g.*prod( (soslp(:,1:3)*omega_vec)./(soslp(:,4:6)*omega_vec),1);
% 
% alpha = [
%   unwrap(angle(Hlp))
%   unwrap(angle(Hhp))
%   unwrap(angle(Hlp+Hhp))
%   unwrap(angle(Hap))
%   ];
% 
% semilogx( ...
%     omega./2./pi, unwrap(angle(Hlp)), '.', ...
%     omega./2./pi, unwrap(angle(Hhp)), 'r', ...
%     omega./2./pi, unwrap(angle(Hlp+Hhp)), 'g--', ... 
%     omega./2./pi, unwrap(angle(Hap)), 'k--' ... 
% );
% legend('lp', 'hp', 'lp+hp', 'ap');
