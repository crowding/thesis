function response = CenterSurroundModel(params, data)

% Here are things that GB's model doesn't capture that I want to.
% First, there is a switch between motion repulsion and
% assimilation that is signfinicant in the fact that 

% se

% third, there is an adjustment to be made to the strength of the
% local motion signal, as a function of DX> Computing this will
% require motion energy models.

% Fourth, there really is a difference between number and density,
% and hte model needs to reflect that.

% Fifth, there really is a different

% Fourth, Geoff's model is weird in that it assumes carrier and
% envelope motion have the same vairance. In fact if you tkae the
% tact that you are combining local and global signals multiplying
% them together, then you don't know what the input variance in
% each is -- but I think number versus density gives you leverage
% over this problem.

[content, spacing, dx] = deal(d.spacing, d.content, d.dx)

