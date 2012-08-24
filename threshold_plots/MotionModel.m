function [prob,G,L,Msig] = MotionModel(p,s,c,dx)


G = p.dxG*dx;  %global contribution
L = -p.cL*c+p.sL./s; %local contribution

%w = s./(p.wSlope+s);

w = .5;
Mbar =  w.*G+(1-w).*L;

Msig = p.ss./s;

prob = normcdf(Mbar,0,Msig);

%prob = 1./(1+exp(-Mbar./Msig));