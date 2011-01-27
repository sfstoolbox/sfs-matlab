% simulate broadband signal
% S.Spors / 20.2.2008

%f2=linspace(20,20000,1000);
f2=logspace(log10(20),log10(20000),200);
%dist=0.05;

P2=zeros(length(f2),1);

for m=1:length(f2)
    f=f2(m);
    main;
    P2(m,:)=P;
    m
end


%plot_fresponse;


if(1)
% plot frequency response
figure
semilogx(f2,db(abs(P2)./max(abs(P2(:)))));


axis([20 20000 -30 5]);
%legend('WFS','SDM y_{ps}=1 m','SDM y_{ps}=0.1 m','SDM y_{ps}=0.01 m','Location','SouthEast');
xlabel('frequency -> Hz');
ylabel('normalized magnitude -> dB');
grid on;
end



if(0)
% plot energy
figure
semilogx(f2,db(cumsum(abs(P2).^2)));

axis([20 20000 0 100]);
xlabel('frequency -> Hz');
ylabel('energy -> dB');
grid on
end