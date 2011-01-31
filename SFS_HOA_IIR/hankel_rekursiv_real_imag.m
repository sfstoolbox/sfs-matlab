%==========================================================================
% Berechnung der Hankelfuntkion rekursiv: 
%       j_(n+1)(kr) = 2n+1/kr * j_(n)(kr) - j_(n-1)(kr)  
%==========================================================================

%Hilfsfunktion für Parameterangabe
function [x,config] = hankel_rekursiv_real_imag(config)
q = linspace(1,config.u,config.u);
%f = linspace(1,20000,config.u); %Frequenzvektor mit Auswertungstellen der Hankelfunktion
n = 1; %Beginne bei Ordnung 2 
m = 28; %Berechnung bis zur Ordnung m
a = real(sphbesselh(0,2,config.k2.*config.R)); %Startwert 1: Hankelfunktion 0ter Ordnung
b = real(sphbesselh(1,2,config.k2.*config.R)); %Startwert 2: Hankelfunktion 1ter Ordnung
c = imag(sphbesselh(0,2,config.k2.*config.R));
d = imag(sphbesselh(1,2,config.k2.*config.R));
z = forward(config,q,m,n,a,b); %Berechnung des Realteils der Hankelfunktion
w = forward(config,q,m,n,c,d); %Berechnung des Imaginärteils der Hankelfunktion
x = z+1i*w;
end
function [P] = forward(config,q,m,n,a,b)
P = zeros(m-n+2,length(q)); %Nullmatrix für Hankelfunktionen
P(1,:) = a; %Setze in 1.Zeile H-Fkt. 0ter Ordnung
P(2,:) = b; %Setze in 2.Zeile H-Fkt. 1ter Ordnung
   while n < m
        
         y =  ( 2.*n+1 )./( config.k2.*config.R ).* b-a; %Rekursionsvorschrift
        
         P(n+2,:) = y; %Setze in n-te Zeile n-te H-Fkt.
         n = n+1; % Setze Ordnung der Hankelfunktion 1 höher
         a = b; %Setze Startwert 2 = Startwert 1
         b = y; %Setze Startwert 1 = errechnete Hankelfunktion
         
                
    end
     
end
