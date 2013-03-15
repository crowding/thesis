function var = params2var(params,freeList)
%var = params2var(params,freeList)
%
%Support function for 'fit.m'
%Written by G.M Boynton, Summer of '00

var = cellfun(@(param) params.(param)(:), freeList(:)');

