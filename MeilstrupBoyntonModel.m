function MeilstrupBoyntonModel

% Notes 
% as is hard to fit all without c0 = -.1
% cj is OK
% jb is hard.  No dependency of mu on c
% je is pretty good
% ko is sparse and OK except for s=2.62
% mc is unfittable
% ml is good
% nj is pretty good
% ns is OK
% pbm is pretty good
% sm is sparse
% tl is pretty good
% to is sparse but OK anywa

% Which data sets to fit
expTypes = {'content','spacing','all'};
expType = 'spacing';

subject = 'pbm';

%Initial parameters
switch subject
    case 'as'
        p.sig0 = .065;
        p.sigk = .75;
        p.siga = 3;
        
        p.mua = 2;
        p.mukc = 1.5;
        p.muks = 1.5;
        p.mu0 = -.1;
        
    case 'jb'
        p.sig0 = .1;
        p.sigk = .75;
        p.siga = .5;
        
        p.mua = 1;
        p.mukc = .5;
        p.muks = 0;
        p.mu0 = -.2;
    otherwise
        % These are some pretty good initial parameters for most subjects
%         p.sig0 = .065;
%         p.sigk = 2;
%         p.siga = 11;
%         
%         p.mua = 3.2;
%         p.mukc = 2;
%         p.muks = 7.7;
%         p.mu0 = 0;
        
        p.siga = 3;
        p.sigk = .75;
        p.sig0 = .75;
        
        p.mua = 8;
        p.muks = 1;
        p.mukc = 1;
        p.mu0 = 0;
         
end
freeParams = {'sig0','sigk','siga','mua','mukc','muks'};

%freeParams = {'sigk','siga','mua','mukc','muks'};


% Load the data if it isn't already loaded ('see UnpackPeter.m')
if ~exist('Sdata','var')
    Sdata = load('data');
end

% Copy the approprate structure to 'data'
data = Sdata.([expType 'Data']);

%% Pull out the relevant trials and fit the model to the data
% List of all subjects
subjectList = unique(data.subject);
nSubjects = length(subjectList);

% Find trials for the current subject
id = strcmp(data.subject,subject);

% Pull out all conditions and responses for all trials for this
% subject
sData = structfun(@(x)x(id), data, 'UniformOutput', 0);
spacing = sData.target_spacing;
content = sData.folded_direction_content;
dx = sData.folded_displacement;
response = 1 - sData.folded_response_with_carrier; %note the 1-response.

%Fit the model (uncomment the next line to see initial parameter predictions)
p = fit(@fitMotionModel, p, freeParams, spacing, content, dx, response);
p;

%% Plot each psychometeric function and model prediction
close all

%Find the list of parameters used in the experiments using 'unique'
eccentrities = unique(data.eccentricity(id));  %not used, always one eccentricity (6.667 deg)
[sList, ~, sIndex] = unique(spacing);   
[dxList, ~, dxIndex] = unique(dx);
[cList, ~, cIndex] = unique(content);

colList = cool(length(cList));  %colormap for plotting

%calculate number of samples and probability correct for each
%condition
nc = accumarray([sIndex, cIndex, dxIndex], response);
ncc = accumarray([sIndex, cIndex, dxIndex], 1-response);
n = nc + ncc;
pc = nc ./ n; %this will be NaN

%Loop through each spacing and plot the psychometric functions

%Obtain predictions from the model
model_pred = nan(size(n));
for sNum = 1:length(sList)
    for cNum = 1:length(cList)
        if sum(n(sNum,cNum,:))
            model_pred(sNum, cNum, :) = ...
                MotionModel(p,sList(sNum)*ones(size(dxList)), ...
                            cList(cNum)*ones(size(dxList)), ...
                            dxList);
        end
    end
end

norm_pred = nan(size(n));
%also fit cumulative normals per condition. They go into mu and
%sig, indexed by sNum and cNum.
for sNum = 1:length(sList)
    for cNum = 1:length(cList)
        %Initial parameters for the cumulative normal
        [prob,predMu,predSig] = MotionModel(p,sList(sNum),cList(cNum),0);
        pNorm.sig =  predSig;       
        pNorm.mu = predMu;
        pNorm.shutup = 'yes'; 

        %Find all trials for this s and c (and all dx)
        id = spacing == sList(sNum) & content == cList(cNum);
        if sum(id)
            %Fit the cumulative normal (with slope fixed by model fit?)
            pNormBest = fit(@fitCumNormal,pNorm,{'mu'},dx(id),response(id));
            %Save the best-fitting parameters
            mu(sNum,cNum) = pNormBest.mu;
            sig(sNum,cNum) = pNormBest.sig;
        end
        norm_pred(sNum, cNum, :) = normcdf(dxList,mu(sNum,cNum),sig(sNum,cNum));
    end
end

%then plot the whole mess
for sNum = 1:length(sList)
    clear h
    figure(sNum)
    clf
    hold on
    title(sprintf('subject %s, s = %5.2f',subject,sList(sNum)));
    % loop through the c's, generating separate psychometric
    % functions
    for cNum = 1:length(cList)
        % loop through the dx's
        for dxNum = 1:length(dxList)
            %plot that point, scaled in size if there are any trials here
            if n(sNum,cNum,dxNum)  
                plot(dxList(dxNum),pc(sNum,cNum,dxNum),'o',...
                    'MarkerSize',n(sNum,cNum,dxNum)/20+5,...
                     'MarkerFaceColor',colList(cNum,:));
            end
        end
        %plot a dashed line through the points
        h(cNum) = plot(dxList,squeeze(pc(sNum,cNum,:)),':','Color',colList(cNum,:));

        if sum(n(sNum, cNum, :))
            %plot the model prediction as a solid line
            plot(dxList,squeeze(model_pred(sNum, cNum, :)),'-', ...
                 'Color',colList(cNum,:));
            %Plot the best fitting cumulative normal
            plot(dxList,squeeze(norm_pred(sNum, cNum, :)),'-.', ...
                 'Color',colList(cNum,:));
        end
    end
    xlabel('dx');
    ylabel('p(clockwise)');
    legend(h,num2str(cList),'Location','SouthEast');
end




%% Plot the model parameters
figure(14);
clf

subplot(2,1,1) %mu as a function of s and c
hold on
clear h
for cNum =1:length(cList)
    %plot mu from the best fitting Cumulative Normal 
    h(cNum) = plot(sList,mu(:,cNum),'o','Color',colList(cNum,:));
    %get and plot mu from the model fit 
    [prob,predMu,predSig] = MotionModel(p,sList,cList(cNum),0);
    plot(sList,predMu,'-','Color',colList(cNum,:));
end
xlabel('s');
ylabel('mu');
legend(h,num2str(cList));

subplot(2,1,2); %sig as a function of s and c
hold on
clear h
for cNum =1:length(cList)
    %Plot sig for the best fitting cumulative normal
    id = sig(:,cNum) >0;
    h(cNum) = plot(sList(id),sig(id,cNum),'o','Color',colList(cNum,:));
    %get and plot sig from the model fit
    [prob,predMu,predSig] = MotionModel(p,sList,cList(cNum),0);
    plot(sList,predSig,'-','Color',colList(cNum,:));
end
xlabel('s');
ylabel('sig');
legend(h,num2str(cList));

%Tile the figures
switch(expType)
  case 'content'
      tile(2,2)
  case 'spacing'
      tile(3,4)
  case 'all'
    tile(3,4)
end
%% plot the underlying model
[svals, cvals, dvals] = ndgrid(.5:.5:21, 0:.05:1, -1:.02:1);
[a1, b1, c1] = MotionModel(p,svals,cvals,dvals);
[probM,muM, sigM] = MotionModel(p,svals,cvals,dvals);

% svals=.5:.5:21;  
% cvals=0:.05:1;
% dvals=-1:.02:1;  
% for s=1:length(svals)
% for c=1:length(cvals)
%     for d=1:length(dvals)
% 
%           [a1, b1, c1] = MotionModel(p,svals(s),cvals(c),dx);
%         [probM(c,s, d),muM(c,s, d), sigM(c, s, d)] = MotionModel(p,svals(s),cvals(c),dvals(d));
%     end
% end
% end

figure(100)
subplot(1, 3, 1)
imagesc(probM(:, :, round(length(dvals)/2)));
colormap(gray)

xlabel('direction content')
ylabel('spacing');
title('probability saying clockwise')

subplot(1, 3, 2)
imagesc(muM(:, :, round(length(dvals)/2)));
colormap(gray)
xlabel('direction content')
ylabel('spacing');
title('relative strength of carrier motion')

subplot(1, 3, 3)
imagesc(sigM(:, :, round(length(dvals)/2)));
colormap(gray)
xlabel('direction content')
ylabel('spacing');
title('standard deviation')

save modelResults.mat