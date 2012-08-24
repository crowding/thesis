%PeterMotionModel.m

p.dxG =4;  %influence of dx on global response G   
p.ss = 1; %influence of s on sensitivity
p.cL = .5;  %influence of c on local response L
p.sL = 3;  %influence of s on local response



sList = [2.09,2.62,3.49,4.65,5.98,8.38,10.47,20.94];

dxList = linspace(-.5,.1,51);
[dx,s] = meshgrid(dxList,sList);

c = .15*ones(size(dx));

prob = MotionModel(p,s,c,dx);

%equivalently:
%M= normcdf(G,L,Msig);

figure(1)
clf

plot(dxList,prob);
legend(num2str(sList'),'Location','NorthWest');
xlabel('dx (envelope step size)');
ylabel('P(clockwise)');
set(gca,'YLim',[0,1]);
hold on
plot([0,0],[0,1],'k-');

%% 'Bias'

%'bias' is the log of proportion 'clockwise' for dx = 0.

cList = linspace(0,1,21);
sList = [2.1,7.0,8.4];

[c,s] = meshgrid(cList,sList);
dx = zeros(size(c));

[prob,G,L,Msig] = MotionModel(p,s,c,dx);



figure(2)
clf
plot(cList,prob)
xlabel('c (directional content)');
legend(num2str(sList'));
ylabel('bias ( P(clockwise|dx = 0) )');

%% Sensitivity

%'sensitivity' is the slope of the psychometric function somewhere

figure(3)
clf
plot(cList,Msig)
xlabel('c (directional content)');
legend(num2str(sList'));
ylabel('Sensitivity');




