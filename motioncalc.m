function trial = motioncalc(trial, mkfig)

r = trial.trial_extra_r;
nt = trial.trial_extra_nTargets;

if ~exist('mkfig', 'var')
    mkfig  = 0;
end

interval = 1/120; 
% I really should pull this from the data but I always
% had 120. If you ask about effects because of randomized onset times on
% the order of less than 1/120 second I'll get stabbity. Although I will
% have to deal with trial by trial effects on the occluded trials.

CP = CauchyPatch...
    ( 'size', [trial.trial_extra_wavelengthScalar*r trial.trial_extra_widthScalar*r trial.trial_extra_durationScalar*trial.trial_extra_dt]...
    , 'order', trial.trial_motion_process_order ...
    , 'velocity', trial.trial_extra_wavelengthScalar*r*trial.trial_extra_tf ...
    , 'phase', 0 ... %a test...
    );

motion = ApparentMotion...
    ( 'primitive', CP ...
    , 'dx', trial.trial_extra_globalVScalar * trial.trial_extra_dt * r ...
    , 'dt', trial.trial_extra_dt ...
    , 'n', trial.trial_motion_process_n ...
    , 'center', [trial.trial_extra_phase * r 0 0] ...
    );

%then replicate the apparent motion around the circle...
stim = ApparentMotion ...
    ( 'primitive', motion ...
    , 'dx', 2*pi*r / trial.trial_extra_nTargets ...
    , 'dt', 0 ...
    , 'n', nt);

%now let's deal with the motion condition (congruent, incongruent, etc.)
switch trial.motionCondition{1}
    case 'congruent'
        contrast = 1/sqrt(2);
    case 'incongruent'
        contrast = 1/sqrt(2);
        stim.primitive.primitive.velocity = -stim.primitive.primitive.velocity;
    case 'ambivalent'
        stim2 = stim;
        contrast = 0.5;
        stim2.primitive.primitive.velocity = -stim2.primitive.primitive.velocity;
    otherwise
        contrast = NaN;
end

%now let's render a picture. We render it for an extent of 2*pi*r in the
%X-direction, and only at 0 in the y-direction.
[x, y, t] = extent(stim);

%how many points to sample around the circle -- the smallest
%wavelengthScalar is 0.05, which would mean the nyquist sampling limit is
%at 4*pi/0.05 = 250 samples. So I'll take 1024 samples to be relatively
%safe and nice to an FFT algorithm. This is a more than pixels on
%the screen in the smallest eccentricities, but I don't think
%there are any resolution problems with the display.
np = 1024;
nt = 256;

%generate the grid. The y-coordinate is just sampled at 0.
xs = x(1):(2*pi*r/np):max(x(2), x(1)+(2*pi*r));
ys = 0;
ts = t(1):interval:max(t(end), t(1)+interval*nt-1);

%evaluate the motion.
z = evaluate(stim, xs, ys, ts);
if isequal('ambivalent', trial.motionCondition{1})
    z = z + evaluate(stim2, xs, ys, ts);
end
z = z * contrast;
%
%Now wrap things around the circle. we should end up with np points
%exactly.
if size(z,2) < np
    z(1, np, 1) = 0;
end
while size(z, 2) > np
    %this proably breaks on wrapping more than once but that never happens
    %in the experiment
    z(:,(np+1):end,:) = z(:,1:size(z,2) - np,:) + z(:,(np+1):end,:);
    z(:, 1:(size(z,2) - 1024), :) = [];
end

% deal with the occlusion condition
xxs = wrap(xs, 0, 2*pi*r);
switch(trial.visibilityCondition{1})
    case 'left'
        vis = (xxs < (4*pi/12 * r)) | (xxs > (20*pi/12 * r));
        occlFactor = 3;
    case 'right'
        %visible ....
        vis = (xxs > (8*pi/12 * r)) & (xxs < (16*pi/12 * r));
        occlFactor = 3;
    case 'full'
        %do nothing
        vis = ones(size(xxs));
        occlFactor = 1;
    otherwise
        trial.visibilityCondition{1}
end
z(:,~vis,:) = -0.05; %the occlusion contrast...

%ditch the useless first dimension...
z = shiftdim(z, 1);
xs = xs(1:size(z,1));



%now let's filter the image.

%we'l use four matched filters according to the four values of
%wavelengthMultiplier that were in the experiment...
%(0.0500 0.0750 0.1125 0.1688)

%these are the filter sizes chosen by Banks, Sekuler and Anderson 1991,
%which I will follow. Note there is some aggrivation in how to match a
%Cauchy filter to a Gabor (matching the bandwidth leads to a Cauchy filter
%that looks too small, but on the other hand Banks et al. make a very
%generous estimate of the region of simple summation. Their intent was to
%have a stimulus size that was for sure large enough to fill all the simple
%summation in the channels, so it's an overestimate of RF size. Right now,
%I have a WAG that the two effects cancel each other out, but I'd like to
%know what kind of variation in receptive field size with spatial frequency
%has been seen physiologically (or with a psychophisical experiment more
%attuned to an unbiased estimate) Where's the physiological literature?

%Real and imaginary parts of 4 filters. see filtersize.m for how I choose
%the bandwidth of each filter.
wmults = [0.0500 0.0750 0.1125 0.1688];
[filtersizes, filtercycles] = arrayfun(@filtersize, ...
                                       repmat(r, size(wmults)), 1./(r.*wmults));

phases = [0 pi/2];
[wmults, phases] = ndgrid(wmults, phases);
spfilt = arrayfun...
    (@(wmult, n, phase) ...
        CauchyPatch( 'size', [wmult*r 1 1]...
                   , 'order', n...
                   , 'phase', phase)...
    , wmults, filtersizes' * [1 1], phases, 'UniformOutput', 0);


%now construct the temporal filters. I'll use the same ones as in Kiani et
%al 2008.

tfilt{1} = @(t) (60*t).^3 .* exp(-60*t).*(1/6   - (60*t).^2/(120) );
tfilt{2} = @(t) (60*t).^5 .* exp(-60*t).*(1/120 - (60*t).^2/(5040));
%now, Kiani et al chose these filters to be consistent with MT temporal
%responses (citing a Movshon ARVO abstract for this; useless) and "centered
%on the speed of stimuli" for their display. Most of my stimuli are
%centered on 10 Hz, and have no reason to believe that I need to use more
%than one temporal filter. e.g. in the TF experiment, but I'll take a broad
%bandwidth. I don't think that centering on the speed makes much of a
%difference.

%now compute each filter on the same size grid as the stimulus
xc = xs - mean(xs); %center the filter
spfilt = cellfun(@(filt) evaluate(filt, xc, 0, 0), spfilt, 'UniformOutput', 0);
tc = ts - min(ts); %temporal filter only coherently defined for t > 0
tfilt = cellfun(@(filt) filt(tc), tfilt, 'UniformOutput', 0);

%normalize all the spatial filters for unit energy. Should the outputs
%undergo probability summation afterwards, perhaps? Adjusted for the
%typical contrast thresholds? Best to move back to R for crunching after
%we're done convoluting.
spfilt = cellfun(@(x) x / norm(x), spfilt, 'UniformOutput', 0);

%TODO also write down a shpiel to narrow down this argument that we keep having,
%so that it can go in the paper and we don't need to keep having it -- the
%justification for looking at physiological-like passbands.

%now that we have the stimulus and all the filters, convolve!!!! via
%fourier transforms, and seperably. basically, take the fourier transform
%of everything, then the cartesian product of 'em all, then sum of squared
%magnitudes for motion energy.

%now, better to work in complex numbers for the separable filters part.
%temporal filter becomes 'forwards' and 'backwards' temporal filter.
spfilt = cellfun(@(r, i) complex(r, i), spfilt(:,1), spfilt(:,2), 'UniformOutput', 0);
tfilt = {complex(tfilt{1}, tfilt{2}), complex(tfilt{1}, -tfilt{2})};

%this gives eight filters (4 x 2).

%filter in the fourier domain, with a cartesian elementwise product....
zf = fft2(z);
spfiltf = cellfun(@fft, spfilt, 'UniformOutput', 0);
spfiltf = cellfun(@transpose, spfiltf, 'UniformOutput', 0);
tfiltf = cellfun(@fft, tfilt, 'UniformOutput', 0);
bsxtimes = @(a, b) bsxfun(@times, a, b);
filtered = cartcellfun(bsxtimes, cartcellfun(bsxtimes, {zf}, spfiltf), tfiltf);


if (mkfig)
    %TODO note from progress meeting: let's put units on the fourier transform of
    %the stimulus and just look at what it does (in a reasonable region) as the
    %stimulus properties are varies. maybe make a demo with knobs to turn.

    %show the stimulus (in natural units) and its fourier transform (in
    %sensible units.)
    subplot(3, 1, 1);
    image(xs, ts, repmat(z'/2+0.5, [1 1 3])); xlabel('degrees (around circle)'); ylabel('sec');
    %set the aspect ratio so that the nominal 100ms delta-x and 0.05*r
    %delta-t are equal.
    daspect([0.05*r 0.1 1]);
    
    subplot(3, 2, 3);
    
    %put units of Kz and cyc/deg on this, and maybe plot the bandwidth of
    %each filter? Note that there's a lot or fencepost errors you can make
    %calculating FFT frequencies. Draw it out on paper.
    zf = fft2(z');
    
    %the length covered by the x and t-sampling (cyclic!)
    lx = numel(xs)*(xs(end)-xs(1))/ (numel(xs)-1);
    lt = numel(ts)*(ts(end)-ts(1)) / (numel(ts)-1);
   
    %what physical frequencies do the fft slots correspond to
    xfs = (0:numel(xs)-1) .* 1/lx; %in cyc/degree
    tfs = (0:numel(ts)-1) .* 1/lt; %in cyc/degree

    %flop the top end onto more sensible negative frequencies.
    tfs = wrap(tfs, -numel(ts)/lt/2, numel(ts)/lt);
    xfs = wrap(xfs, -numel(xs)/lx/2, numel(xs)/lx);    
    
    zf = fftshift(zf);
    tfs = fftshift(tfs);
    xfs = fftshift(xfs);
    
    %the range of frequencies we are interested in: -20 to 20 cpd and -40 to
    %40 Hz, say..
    
    keeptfs = abs(tfs) < 20;
    keepxfs = abs(xfs) < 10;
    zf = zf(keeptfs, keepxfs);
    ktfs = tfs(keeptfs);
    kxfs = xfs(keepxfs);
    
    %let's chop out what's in the stimulus, below 20 cps and below 5
    %cyc/degree, on a decibel scale.
    %note that I downscale the FFT image because MATLAB doesn't do it right
    %when displaying, eading to artifacts.
    
%    imagesc...
%        ( linspace(kxfs(1), kxfs(end), 128)...
%        , linspace(ktfs(1), ktfs(end), 128)...
%        , imresize(20*log10(abs(zf)), [128 128])...
%        , [-40 30]);
    imagesc...
        ( linspace(kxfs(1), kxfs(end), 128)...
        , linspace(ktfs(1), ktfs(end), 128)...
        , imresize(abs(zf), [128 128])...
        );
    axis square;
    colormap(gca, hot(256));
    xlim([-10 10]);
    ylim([-20 20]);
    %cba = colorbar;
    %ylabel(cba, 'dB');
    
    %now also show the FWHM of the various filters on top.
    hold on;
    
    %ugh, feed in linespecs....
    colors = reshape({'m','g', 'c', 'w'}, size(spfiltf));
    lines = reshape({'-', ':'}, size(tfiltf));
    
    cartcellfun(@plotfwhm ...
        , cellfun(@(a,b){a,b},spfiltf,colors,'UniformOutput', 0)...
        , cellfun(@(a,b){a,b},tfiltf,lines,'UniformOutput', 0)...
        );
    hold off;
    xlabel('spatial freq. (cyc/deg)'); ylabel('temporal freq(Hz)');
    set(gca, 'YDir', 'normal');
    %todo: show the FWHM of the various filters on top (contour plot?)
end

function out = plotfwhm(s,t)
    sss = s{1};
    ttt = t{1};
    ls = [s{2} t{2}];
    spec = abs(fftshift(sss * ttt));
    contour(xfs, tfs, spec'/max(max(spec)), [0.5 0.5], ls);
    out = [];
end

%sum up the energy in each channel(filter response squared.)
energies = cellfun(@(x) sum(sum(real(x.*conj(x)))), filtered);
energies = reshape(energies, [numel(spfiltf) numel(tfiltf)]);

%normalize spectra per number of elements?
nelem = trial.trial_motion_process_n * trial.trial_extra_nTargets;

%this should just be the same up to a factor (norms are preserved under orthonormal change of basis), but I'll double check
%energiesf = cellfun(@(x) sum(sum(real(x.*conj(x)))), filtereds);
%energiesf = reshape(energiesf, [numel(spfiltf) numel(tfiltf)]);

%let's see each of the filtered responses...
if (mkfig)
    subplot(3,2,4);
    barh(1./(r*wmults(:,1)), energies(:,1) / nelem / trial.trial_extra_wavelengthScalar / trial.trial_extra_durationScalar / trial.trial_extra_dt * occlFactor);
    hold on;
    barh(1./(r*wmults(:,2)), -energies(:,2) / nelem / trial.trial_extra_wavelengthScalar / trial.trial_extra_durationScalar / trial.trial_extra_dt * occlFactor);
    hold off;
    title('leftward and rightward energy');
    xlabel('energy');
    ylabel('filter spatial frequency (cyc/s)');
    xlim([-5*10^9, 5*10^9])
end


%record the outputs from each filter in the struct, and we're done for now.
%in a singleton cell (as it'll be converted back into an array of structs)
trial.analysis_wavelengths = {r.*wmults(:,1)'};
trial.filter_sizes = {filtersizes};
trial.filter_cycles = {filtercycles};
trial.energy_global = {energies(:,1)'};
trial.energy_nonglobal = {energies(:,2)'};

%motion energy demands we sum in quadrature...
%now after we've done this, sum up to find the total energy (would it be
%useful to see the energy space/time separated? Ha, not gonna think about
%that.

%split into leftward and rightward filters...

end