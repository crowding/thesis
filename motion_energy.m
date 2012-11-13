function motion_energy(infile, outfile)

% create space-time example stimulus plots of all unique stimuli.
% calculate the motion energy of the given stimuli in each band.  The
% inpus file is the list of all trial parameters.  The output file,
% however, only tracks unique trials (up to a rotation and
% reflection.)

if exist('infile', 'var')
    infile = 'data.mat';
end

S = load(infile);

stimuli = soa2aos(stimuli);

%now, to a motion-energy calculation for each stimulus.
stimuli = arrayfun(@motioncalc, stimuli);

%and the output, back to R
stimuli = aos2soa(stimuli);

save(outfile, 'stimuli');

end