function[pw_b, pw_m] = HOA_driving_signal_broadband_mono_pw(config)
%==========================================================================

% calculation of driving signal broadband plane wave

L=length( config.x0 );
F=zeros( length( config.k2 ), L );  %Matrix for coefficients and LS

al=atan2( config.x0( 2,: ), config.x0( 1,: ) );

% get spatial band limitation 
if mod( L,2 )
   N21 =  floor( L/2 );
   N22 = N21;
else
   N21 = L/2-1;
   N22 = L/2;
end


for l=0:L-1 % loop over all loudspeakers   
    
  disp([l L-1]); %display computation time of loop
    
    for m = -N21 : N22
       
        F(:, l+1) = F(:, l+1) + (2./config.R)* -1i^(abs( m ) )./ (-1i*config.k2) * 1./ ...
            sphbesselh(abs( m ), 2, config.k2.*config.R) .* exp(1i*m* (config.al_pw-al( l+1 ) ) );
   
    end
    
end
%=========================================================================
% gd = zeros(length(config.k2),L);
% 
% for n = 0:1:L-1
%    
%     disp([n L-1])
%         
%     grp_delay = grpdelay(F(:,n+1),1,config.k2,config.fs);
%     gd(:,n+1) = gd(:,n+1) + grp_delay;
%     
% end

%==========================================================================
%calculation of driving signal for pw

L=length( config.x0 );
E=zeros( 1,L );


for l=0:L-1 % loop over all loudspeakers   
    
    disp([l L-1]); %diplay computation time of loop
    
    for m = -N21 : N22
    E(l+1) = E(l+1) + (2./config.R)* -1i^(abs( m ))./ (-1i*config.k) * 1./ ...
            sphbesselh(abs( m ), 2, config.k.*config.R) .* exp(1i*m*( config.al_pw-al( l+1 ) ) );
    end
end
%=========================================================================
pw_b = F; % Driving function broadband
%pw_m = 1;
pw_m = E; % Driving function monofrequent
end
%=========================================================================