function[ps_b_r,ps_m_r] = HOA_driving_signal_broadband_mono_ps_rek(config)
%=========================================================================
z = hankel_rekursiv(config); %rekursiv calculation of Hankelfunctions with r = R
h2 = hankel_rekursiv_ps(config); %rekursiv calculation of Hankelfunctions with r = rs_rel
%=========================================================================
%Calculation driving signal broadband point source
L=length(config.x0);
E=zeros(length(config.k2), L);  %Matrix for coefficients and LS

% get spatial band limitation 
if mod(L,2)
   N21 =  floor(L/2);
   N22 = N21;
else
   N21 = L/2-1;
   N22 = L/2;
end

% Position of sound source (relative to center of LS-array)
xs_rel = ( config.xps(1) - config.xref(1) );
ys_rel = ( config.xps(2) - config.xref(2) );


rs_rel = sqrt(  xs_rel.^2 + ys_rel.^2 ); %  get distance from point source
alphas_rel = atan2( ys_rel, xs_rel );   % get angle from position of point source

h2 = h2.';
z = z.';

for l=0:L-1 % loop over all loudspeakers   
    
    alpha_0 = 2*pi/L * l;   % get angle for every LS in relation to given point
    
    disp([l L-1]);  % show time of computation for coeffcients
    
    for m = -N21 : N22
     
        E(:, l+1) = E(:, l+1) +  1/(2*pi*config.R) .* h2(:,abs(m)+1) ./ ...
            z(:,abs(m)+1) .* exp( 1i * m * (alpha_0 - alphas_rel) );
    
    end
    
end
%%=========================================================================

% % calculation of driving signal for ps monofreqent
L=length(config.x0);
D=zeros(1, L);  %Matrix for coefficients and LS


for l=0:L-1 % loop over all loudspeakers   
    
    alpha_0 = 2*pi/L * l;   % get angle for every LS in relation to given point
    
    disp([l L-1]);  % show time of computation for coeffcients
    
    for m = -N21 : N22
     
        D(l+1) = D(l+1) +  1/(2*pi*config.R) .* sphbesselh(abs(m),2, config.k.*rs_rel) ./ ...
            sphbesselh(abs(m),2, config.k.*config.R) .* exp( 1i * m * (alpha_0 - alphas_rel) );
    
    end
%     

%==========================================================================

%==========================================================================
ps_b_r = E; % Driving function broadband
% ps_m_r = 1;
ps_m_r = D; % Driving function monofrequent
%==========================================================================
end
