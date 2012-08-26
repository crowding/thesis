function [params,err] = fit(fun,params,freeList,varargin)
%[params,err] = fit(@fun,params,freeList,var1,var2,var3,...)
%
%Helpful interface to matlab's 'fmins' function.
%
%INPUTS
% @fun     :  function be optimized.  Has form fun(params,var1,var2,...)
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
%disp(sprintf('Fitting "%s" with %d free parameters.',fun,length(var)));
var = fminsearch(@fitFun,var,options,fun,params,freeList,varargin);

%get final parameters
params=  var2params(var,params,freeList);

%evaluate the function
err = fun(params, varargin{:});

return

%debug

clear params








