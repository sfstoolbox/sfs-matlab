% simulate wavefield reproduced by a sound reproduction system
% S.Spors, 23.7.2007

%clear all

% simulation parameters

nLS = 72;                  % number of loudspeakers
selectLS = 1;               % apply secondary source selection

al_pw = pi/2;               % incidence angle of plane wave
xps = [0 3];
xfs = [-0.5 0];
xref = [0 0];

% simulation grid
%x = linspace(-2,2,400);
%y = linspace(-2,2,400);
%y = linspace(-0.05,3,1000);

x=0;y=0;



% frequency
%f = 3000;
k = (2*pi*f)/343;


%================================================================================
% get loudspeaker positions and normal vectors

% circular array
R=1.50;
[x0,n0] = LSpos_circ(0,0,R,nLS);
%n0(43)=90/180*pi;

% rectangular array
%[x0,n0] = LSpos_box(0,0,3,nLS);

% linear array
%dist=0.15;
%[x0,n0] = LSpos_linear(0,0,(nLS-1)*dist,nLS);


%================================================================================
% compute loudspeaker driving signals

% 3D WFS
%[D] = driving_signal_pw(al_pw,k,x0,n0,selectLS);
%[D] = driving_signal_ps(xps,k,x0,n0,selectLS);
%[D] = driving_signal_fs(xfs,k,x0,n0,selectLS,-pi/2);

% 2.5D WFS
%[D] = WFS25D_driving_signal_pw(al_pw,k,x0,n0,xref,selectLS);       %D(29)=0;
%[D] = WFS25D_driving_signal_ps(xps,k,x0,n0,xref,selectLS);
%[D] = WFS25D_driving_signal_ps_TUD(xps,k,x0,n0,xref,selectLS);
%[D] = WFS25D_driving_signal_ps_AES(xps,k,x0,n0,xref,selectLS);
[D] = WFS25D_driving_signal_fs(xfs,k,x0,n0,xref,selectLS,0);

% 2.5D local WFS
%[x00,n00] = LSpos_circ(0,-0.5,0.2,nLS);
%[D] = WFS25Dlocal_driving_signal_pw(al_pw,k,x0,n0,x00,n00,selectLS);


% 2D WFS
%[D] = driving_signal_pw(al_pw,k,x0,n0,selectLS);
%[D] = driving_signal_ls(xps,k,x0,n0,selectLS);
%[D] = driving_signal_fls(xfs,k,x0,n0,selectLS,pi);
%[D] = driving_signal_ls_OS(xps,k,x0,n0,selectLS,dist);

% oversampled circular 2D WFS
if(0)
   osf=10;
   
   [x0i,n0i] = LSpos_circ(0,0,R,osf*nLS);
   [Di] = driving_signal_pw(al_pw,k,x0i,n0i,selectLS);
   
   D = resample(Di,1,osf);
   D(30:end)=0;
    
end

% 2.5D HOA
%[D] = HOA25D_driving_signal_pw(al_pw,k,x0,R);

% 2D HOA
%[D] = HOA_driving_signal_pw(al_pw,k,x0,R);

% amplitude panning Ambisonics
%[D] = APA_driving_signal_pw(al_pw,k,x0,R);



% 2.5D SDM
%[D] = SDM25D_driving_signal_pw(al_pw,k,x0,n0,xref,selectLS); 
%[D] = SDM25D_driving_signal_ps_approx(xps,k,x0,n0,xref,selectLS);



% rectangular SFR
%[D] = DFT_driving_signal_pw(al_pw,k,x0,n0);


%================================================================================
% apply tapering window

if(1)
    idx=find(abs(D) > 0);
    L=length(idx)*0.10;
    win = hanningwin(L,L,length(idx));
    D(idx)=D(idx).*win';
end


%================================================================================
% compute reproduced wave field

P = zeros(length(x),length(y));

idx=find(abs(D)~=0);
for n=idx
   P = P + D(n) .* point_source(x,y,x0(1,n),x0(2,n),k);
   %P = P + D(n) .* line_source(x,y,x0(1,n),x0(2,n),k);
end



%================================================================================
% plot reproduced wave field

figure; plot_wavefield;