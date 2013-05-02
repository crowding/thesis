function err = fitFun(var,fun,params,freeList,origVarargin)
%err = fitFun(var,@fun,params,freeList,origVarargin)
  
%stick values of var into params

params = var2params(var,params,freeList);

%evaluate the function
err = fun(params, origVarargin{:});

return



