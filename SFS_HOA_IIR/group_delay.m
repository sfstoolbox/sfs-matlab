function [delay] = group_delay(config, E)

%calculate delay line for the desired filter with group delay of the driving function 

var1 =size(E); % get size of matrix E
var2 = zeros(var1(1)-1,config.nLS); %create matrix one line shorter then E due to diff() function
delay = zeros(1,56);    % create vector for delay line

for i = 1:config.nLS % calculate group delay for every LS
    
    phase = unwrap(angle(E(:,i))); % get desired phase 
    %d(phase)/d(omega) =  -(phase(omega+d(omega))-phase(omega))/(d(omega))
    var2(:,i) = -(diff(phase)*config.u/(2*pi*config.fs)); % calculate group delay of phase:
    delay(1,i) = 100000*var2(1,i); % expand delay with factor x due to buffer shifting is only possible with integer values
    
end

delay = round(delay); % round delay to an integer
 
end