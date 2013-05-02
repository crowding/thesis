function params = var2params(var,params,freeList)
%params = var2params(var,params,freeList)
%
%Support function for 'fit.m'
%Written by G.M Boynton, Summer of '00

count = 1;
for i=1:length(freeList)
    len = length(params.(freeList{i}));
    params.(freeList{i}) = var(count:count+len-1);
    count = count+len;
end

