function [prob,mu,sig,vars] = MotionModel(p,d)
%[prob,mu,sig] = MotionModel(p,d)
%
% Inputs:
% p structure containing model parameters
% d a struct with fields:
% 'spacing'
% 'contentt'
% 'dx'
%
% ---- notes from PBM ___ The following are the model
% parameters. These appear to be controlling the shift of two
% Gaussians which are being multiplied together. So we have
% coefficients for the shift of mean and standard deviation for both.
%
% p.mua    % weighting factor for
% p.siga
% p.mukc
% p.sigkc
% p.muks
% p.mu0
%
% One think that the motivation of this model assumes is that the
% motion signal is represented as a continuous variable, The other
% is that the
%
% An interesting question is, in how much of the fits do the Gaussians
% even overlap at all? If the d-prime between the two signals is
% particularly large, then we are really relying the (e^-(x^2)) decay
% of the tail of a Gaussian, rather than this shape. If it is doing
% this, it would explain why there are other exponentials in the model
% which make no sense to me.
%
% For example: an exponential on motion content makes no sense; motion
% content is somewhat analagous to "coherence" in a random dot
% stimulus which you usually take the *log* of to come out with
% something that looks like a psychometric function.
%
% Anyhow, this particular model should easily be captured by a
% generalized linear model with binomial error and Gaussian link. That
% would have the advantage of (a) giving you some damn error bars and
% (b) guaranteed (and fast) convergence with IRLS type algorithm
% rather than getting stuck in false minima with fminsearch.
%
% I think I will have to plot what this model makes out of typical
% conditions.
%
% __________________________________________________
%
% Model:
%   Assume that on any trial, there is a global motion signal and local
%   signal that is drawn from two normal distributions _with equal
%   variance_.  A subject responds 'clockwise' when the draw from the global
%   motion signal exceeds the draw from the local motion signal.
%
%   Assume that the mean global signal is proportional to dx.
%
%   The mean of the local motion signal, mu, is a function of c and s. It
%   is a separable function consisting of an exponential that increases
%   with c with a constant p.cck and and exponential that decreases with s
%   with a constant p.ck, scaled by p.a plus an offset p.c0
[s, c, dx] = deal(d.spacing, d.content, d.dx);

mu = -p.mua.*exp(p.mukc*c).*exp(-p.muks*s)+p.mu0;

%   The standard deviation of the two signals decreases with s as an
%   exponential function with constant p.sk, scale p.sa and offsed p.sig0.

sig = p.siga*exp(-p.sigk*s)+p.sig0;

%    This is the probability of responding 'clockwise' on any trial.
prob = normcdf(dx,mu,sig);
