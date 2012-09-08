function params = initialParams(data, splits)

if ~exist('splits', 'var')
    splits = {'subject'};
end

% Notes
% as is hard to fit all without c0 = -.1
% cj is OK
% jb is hard.  No dependency of mu on c
% je is pretty good
% ko is sparse and OK except for s=2.62
% mc is unfittable
% ml is good
% nj is pretty good
% ns is OK
% pbm is pretty good
% sm is sparse
% tl is pretty good
% to is sparse but OK anywa

p = struct();
p.subject = {'as'};
p.sig0 = .065;
p.sigk = .75;
p.siga = 3;

p.mua = 2;
p.mukc = 1.5;
p.muks = 1.5;
p.mu0 = -.1;

param_table = dataset(p);

p = struct();
p.subject = {'jb'};
p.sig0 = .1;
p.sigk = .75;
p.siga = .5;

p.mua = 1;
p.mukc = .5;
p.muks = 0;
p.mu0 = -.2;

param_table = [param_table;dataset(p)];

p = struct();
p.siga = 3;
p.sigk = .75;
p.sig0 = .75;

p.mua = 8;
p.muks = 1;
p.mukc = 1;
p.mu0 = 0;

param_defaults = dataset(p);

%----------------------------------------

conditions = unique(data(:,splits));

[params, ~, ib] = join(param_table, conditions, ...
                       'type', 'inner', 'MergeKeys', true);

%any conditions that weren't matched get default...
%Ugh, doing it thisway was more trouble than it was worth.
unmatched = conditions(~ismember(1:size(conditions, 1), ib), :);
unmatched.dsfargeg = ones(size(unmatched, 1), 1);
param_defaults.dsfargeg = ones(size(param_defaults, 1), 1);

defaulted = join(param_defaults, unmatched, 'type', 'inner', 'MergeKeys', true);
defaulted(:,'dsfargeg') = [];

params = [params; defaulted];

end

