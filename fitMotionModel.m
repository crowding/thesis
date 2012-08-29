function err = fitMotionModel(p,data)
% err = fitMotionModel(p,s,c,dx,response)
% Finds the predicted probabilty clockwise for each trial and calculates
% the -log likelihood that the responses came from this model.

prob = MotionModel(p,data);
prob = prob*.99+.005;
err = -sum(data.response.*log(prob) + (1-data.response).*log(1-prob));