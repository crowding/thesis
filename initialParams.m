function p = initialParams(subject)

switch subject
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
