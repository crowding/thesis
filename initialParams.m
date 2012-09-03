function p = initialParams(data, splits)

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

switch data.subject{1}
    case 'as'
        p.sig0 = .065;
        p.sigk = .75;
        p.siga = 3;

        p.mua = 2;
        p.mukc = 1.5;
        p.muks = 1.5;
        p.mu0 = -.1;

    case 'jb'
        p.sig0 = .1;
        p.sigk = .75;
        p.siga = .5;

        p.mua = 1;
        p.mukc = .5;
        p.muks = 0;
        p.mu0 = -.2;
    otherwise
        % These are some pretty good initial parameters for most subjects
%         p.sig0 = .065;
%         p.sigk = 2;
%         p.siga = 11;
%
%         p.mua = 3.2;
%         p.mukc = 2;
%         p.muks = 7.7;
%         p.mu0 = 0;

        p.siga = 3;
        p.sigk = .75;
        p.sig0 = .75;

        p.mua = 8;
        p.muks = 1;
        p.mukc = 1;
        p.mu0 = 0;
end

