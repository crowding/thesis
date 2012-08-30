function response = CenterSurroundModel(params, data)

% Here is a model more in tune with the center-surround nature of processing.

% I do not believe that the brain is implementing an estimate of
% velocity in this illusion. Rather I would view it as being part
% of a decision process. 

[content, spacing, dx] = deal(d.spacing, d.content, d.dx)