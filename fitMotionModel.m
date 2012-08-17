function err = fitMotionModel(p,s,c,dx,response)
% err = fitMotionModel(p,s,c,dx,response)
% Finds the predicted probabilty clockwise for each trial and calculates
% the -log likelihood that the responses came from this model.

prob = MotionModel(p,s,c,dx);
prob = prob*.99+.005;
err = -sum(response.*log(prob) + (1-response).*log(1-prob));