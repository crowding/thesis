function [params,err] = fit(funName,params,freeList,varargin)
%[params,err] = fit(funName,params,freeList,var1,var2,var3,...)
%
%Helpful interface to matlab's 'fmins' function.
%
%INPUTS
% 'funName':  function to be optimized.  Has form <funName>(params,var1,var2,...)
% params   :  structure of parameter values for fitted function
% freeList :  List of parameter names to be free in fit
% var<n>   :  extra variables to be sent into fitted function
%
%OUTPUTS
% params   :  structure for best fitting parameters
% err      :  error value at minimum
%
%Example:
% params.a = [1.1,1.2,1.3];
% params.b = 2;
% params.c = 3;
% freeList = {'a(1:2)','c'};
%
%[params,err] = fit('myTestFun',params,freeList,fee,fi,fo,fum);
%
%note the use if indices to set a subset of parameters free 'a(1:2)'
%
%Written by gmb, Summer of '00

%turn free parameters in to 'var'


if isfield(params,'options')
  options = params.options;
else
  options = [];
end


if isempty(freeList)
  freeList = fieldnames(params);
end

var = params2var(params,freeList);
%disp(sprintf('Fitting "%s" with %d free parameters.',funName,length(var)));
var = fminsearch('fitFun',var,options,funName,params,freeList,varargin);

%get final parameters
params=  var2params(var,params,freeList);

%evaluate the function

evalStr = sprintf('err = %s(params',funName);
for i=1:length(varargin)
  evalStr= [evalStr,',varargin{',num2str(i),'}'];
end
evalStr = [evalStr,');'];
eval(evalStr);


return

%debug

clear params








