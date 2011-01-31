%==========================================================================
% Berechnung der Hankelfuntkion rekursiv: 
%       j_(n+1)(kr) = 2n+1/kr * j_(n)(kr) - j_(n-1)(kr)  
%==========================================================================

%Hilfsfunktion für Parameterangabe
function [z,config] = hankel_rekursiv_show(config)
q = linspace(1,20,20);
%f = linspace(1,20000,config.u); %Frequenzvektor mit Auswertungstellen der Hankelfunktion
n = 1; %Beginne bei Ordnung 2 
m = 28; %Berechnung bis zur Ordnung m
x = 0:0.001:20;
a = sphbesselh(0,2,x); %Startwert 1: Hankelfunktion 0ter Ordnung
b = sphbesselh(1,2,x); %Startwert 2: Hankelfunktion 1ter Ordnung
z = forward(config,x,m,n,a,b); %Berechnung der Hankelfunktion

end
function [P] = forward(config,q,m,n,a,b)
x=0:1:20000;
P = zeros(m-n+2,length(q)); %Nullmatrix für Hankelfunktionen
P(1,:) = a; %Setze in 1.Zeile H-Fkt. 0ter Ordnung
P(2,:) = b; %Setze in 2.Zeile H-Fkt. 1ter Ordnung
   while n < m
        
         y =  ( 2.*n+1 )./( x ).* b-a; %Rekursionsvorschrift
        
         P(n+2,:) = y; %Setze in n-te Zeile n-te H-Fkt.
         n = n+1; % Setze Ordnung der Hankelfunktion 1 höher
         a = b; %Setze Startwert 2 = Startwert 1
         b = y; %Setze Startwert 1 = errechnete Hankelfunktion
         
                
    end
     
end
