function []=set_angleaxis(waxis)


if(nargin==0)
  waxis = 0;
end


ticks=0:90:360;

switch waxis
 
 case 0 
  set(gca,'XTick',ticks);
  
 case 1
  set(gca,'YTick',ticks);
  
 case 2
  set(gca,'XTick',ticks);
  set(gca,'YTick',ticks);
  
 otherwise
  disp 'no valid option!'
  
end

  
  
  