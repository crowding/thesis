function err = fitCumNormal(p,x,response)
% err = fitMotionModel(p,s,c,dx,response)

prob = normcdf(x,p.mu,p.sig);
prob = prob*.99+.005;
err = -sum(response.*log(prob) + (1-response).*log(1-prob));