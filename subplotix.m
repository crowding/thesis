function h = subplotix(m,n,i,j, varargin)
    %actually get a subplot by rows and columns
    ix = sub2ind([n m], j, i);
    h = subplot(m,n,ix, varargin{:});
end
