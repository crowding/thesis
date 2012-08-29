function [prob,mu,sig] = MotionModel(p,d)
%[prob,mu,sig] = MotionModel(p,data)
%
% Inputs:
% p structure containing model parameters 
% d a struct with fields:
% 'spacing'
% 'contentt'
% 'dx'
% 
%
% Model:
%   Assume that on any trial, there is a global motion signal and local
%   signal that is drawn from two normal distributions with equal
%   variance.  A subject responds 'clockwise' when the draw from the global
%   motion signal exceeds the draw from the local motion signal. 
%
%   Assume that the mean global signal is proportional to dx. 
%
%   The mean of the local motion signal, mu, is a function of c and s. It
%   is a separable function consisting of an exponential that increases
%   with c with a constant p.cck and and exponential that decreases with s
%   with a constant p.ck, scaled by p.a plus an offset p.c0
[s, c, dx] = deal(d.spacing, d.content, d.dx);

mu = -p.mua*exp(p.mukc*c).*exp(-p.muks*s)+p.mu0;

%   The standard deviation of the two signals decreases with s as an
%   exponential function with constant p.sk, scale p.sa and offsed p.sig0.

sig = p.siga*exp(-p.sigk*s)+p.sig0;

%    This is the probability of responding 'clockwise' on any trial.
prob = normcdf(dx,mu,sig);