function [C] = separate(v, n)
%function C = separate(v, n)
%given a vector V of positive integers, returns a cell array C, such
%that for each i 1:n, v(C{i}) == i.
%If N is not specified, n is the maximum of V.
%
%That is, this takes the third output of "unique" and actually
%constructs the grouping array you always want but end up fudging with
%a for loop and a logical mask.
%
%Example: If V = [1 2 3 2 3 2 2 2 1 1 2 1], then
%separate(V) = {[1 9 10 12]';[2 4 6 7 8 11]';[3 5]'}
if (~exist('n', 'var'))
    n = max(v(:));
end

spa = sparse(1:numel(v), v, 1:numel(v), numel(v), n);

C = arrayfun(@(x)nonzeros(spa(:,x)), (1:n)', 'UniformOutput', 0);