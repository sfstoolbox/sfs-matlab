% draws loudspeaker symbols
% S.Spors / 27.7.2007
function [] = draw_loudspeakers(x0,n0,LSactive)

nLS = length(n0);

if(nargin<3)
    LSactive = zeros(1,nLS);
else
    hold on;
end

if(length(LSactive)==1)
    LSactive = LSactive*ones(1,nLS);
end


w = 0.10;   % total size of loudspeaker symbol
h = 0.10;


% vertex coordinates
v1 = [-h/2 -h/2 0 0 -h/2 ; -w/2 w/2 w/2 -w/2 -w/2];

v2 = [0 h/2 h/2 0 ; -w/6 -w/2 w/2 w/6];


% draw loudspeakers
for n=1:nLS

    
    R = [cos(n0(n)) -sin(n0(n));sin(n0(n)) cos(n0(n))];
    
    for k=1:length(v1)
        vr1(:,k) = R * v1(:,k);
    end
    
    for k=1:length(v2)
        vr2(:,k) = R * v2(:,k);
    end
    
    % shift
    v01(1,:) = vr1(1,:) + x0(1,n);
    v01(2,:) = vr1(2,:) + x0(2,n);
    
    v02(1,:) = vr2(1,:) + x0(1,n);
    v02(2,:) = vr2(2,:) + x0(2,n);

    if(LSactive(n)>0)
        fc = [(1-LSactive(n)) (1-LSactive(n)) (1-LSactive(n))];
        
        fill(v01(1,:),v01(2,:),fc);
        fill(v02(1,:),v02(2,:),fc);
    else
        h=line(v01(1,:),v01(2,:));
        set(h,'Color','k');
        h=line(v02(1,:),v02(2,:));
        set(h,'Color','k');
    end
        
    

end
% deg = pi/180;
% w = linspace(1,360,1000);
% si = sin(w*deg);
% co = cos(w*deg);
% plot(co,si)
% if(nargin==3)
%     hold off;
% end

