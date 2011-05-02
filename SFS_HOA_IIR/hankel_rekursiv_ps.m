function [h2, config] = hankel_rekursiv_ps(config)


% Position of sound source (relative to center of LS-array)
xs_rel = ( config.xps(1) - config.xref(1) );
ys_rel = ( config.xps(2) - config.xref(2) );
rs_rel = sqrt(  xs_rel.^2 + ys_rel.^2 ); %  get distance from point source

q = linspace(1,config.u,config.u);
%f = linspace(1,20000,config.u); %Frequenzvektor mit Auswertungstellen der Hankelfunktion
n = 1; %Beginne bei Ordnung 2 
m = 28; %Berechnung bis zur Ordnung m
a = sphbesselh(0,2,config.k2.*rs_rel); %Startwert 1: Hankelfunktion 0ter Ordnung
b = sphbesselh(1,2,config.k2.*rs_rel); %Startwert 2: Hankelfunktion 1ter Ordnung
h2 = forward(config,rs_rel,q,m,n,a,b); %Berechnung der Hankelfunktion

end
function [P] = forward(config,rs_rel,q,m,n,a,b)
P = zeros(m-n+2,length(q)); %Nullmatrix für Hankelfunktionen
P(1,:) = a; %Setze in 1.Zeile H-Fkt. 0ter Ordnung
P(2,:) = b; %Setze in 2.Zeile H-Fkt. 1ter Ordnung
   while n < m
        
         y =  ( 2.*n+1 )./( config.k2.*rs_rel ).* b-a; %Rekursionsvorschrift
        
         P(n+2,:) = y; %Setze in n-te Zeile n-te H-Fkt.
         n = n+1; % Setze Ordnung der Hankelfunktion 1 höher
         a = b; %Setze Startwert 2 = Startwert 1
         b = y; %Setze Startwert 1 = errechnete Hankelfunktion
         
                
    end
     
end