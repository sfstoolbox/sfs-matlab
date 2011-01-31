function[pw_b_r, pw_m_r] = HOA_driving_signal_broadband_mono_pw_rek(config)
%%=========================================================================
%%rekursiv calculation of Hankelfunctions
%%=========================================================================
z = hankel_rekursiv(config); 
% %========================================================================
% % calculation of driving signal pw (broadband)
% %========================================================================
L=length( config.x0 );
F=zeros( length( config.k2 ), L );  %Matrix for coefficients and LS

% get spatial band limitation 
al=atan2( config.x0( 2,: ), config.x0( 1,: ) );

if mod( L,2 )
   N21 =  floor( L/2 );
   N22 = N21;
else
   N21 = L/2-1;
   N22 = L/2;
end

z = z.';



for l=0:L-1 % loop over all loudspeakers   
    
  disp([l L-1]); %display computation time of loop
    
    for m = -N21 : N22
       
        F(:, l+1) = F(:, l+1) + (2./config.R)* -1i^(abs( m ) )./ (-1i*config.k2) * 1./ ...
            z(:,abs(m)+1) .* exp(1i*m* (config.al_pw-al( l+1 ) ) );
   

        
    end
    
end
% %========================================================================
% % calculation of driving signal for pw (monofrequent)
% %========================================================================

 
        
% L=length( config.x0 );
% E=zeros( 1,L );
% 
% for l=0:L-1 % loop over all loudspeakers   
%     
%   disp([l L-1]); %diplay computation time of loop
%     
%     for m = -N21 : N22
%           
%             E(l+1) = E(l+1) + (2./config.R)* -1i^(abs( m ))./ (-1i*config.k) * 1./ ...
%             sphbesselh(abs( m ), 2, config.k.*config.R) .* exp(1i*m*( config.al_pw-al( l+1 ) ) );
%    
%     end
%     
% end

%%=========================================================================
pw_b_r = F; % Driving function broadband
pw_m_r = 1; % choose, if you don't need monofrequent calculation
% pw_m_r = E; % Driving function monofrequent
%%=========================================================================
end