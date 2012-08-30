function MeilstrupBoyntonModel(export)
%if given value of 1, will assign everything to base workspace to
%play with.

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

p = initialParams(subject);

freeParams = {'sig0','sigk','siga','mua','mukc','muks'};

% Load the data if it isn't already loaded ('see UnpackPeter.m')
S = load('data');

% Copy the approprate structure to 'data'
data = S.([expType 'Data']);

%% Pull out the relevant trials and fit the model to the data
% List of all subjects
subjectList = unique(data.subject);
nSubjects = length(subjectList);

% Find trials for the current subject
id = strcmp(data.subject,subject);

% Pull out all conditions and responses for all trials for this
% subject
data = structfun(@(x)x(id), data, 'UniformOutput', 0);
data = rename(data, ...
              'target_spacing', 'spacing', ...
              'folded_direction_content', 'content', ...
              'folded_displacement', 'dx', ...
              'folded_response_with_carrier', 'response');

data.response = 1 - data.response;

%Fit the model (uncomment the next line to see initial parameter predictions)
p = fit(@fitMotionModel, p, freeParams, data);

%% Plot each psychometeric function and model prediction
close all

%Find the list of parameters used in the experiments using 'unique'
eccentrities = unique(data.eccentricity);  %not used, always 6.667 deg
[sList, ~, sIndex] = unique(data.spacing);
[dxList, ~, dxIndex] = unique(data.dx);
[cList, ~, cIndex] = unique(data.content);

colList = cool(length(cList));  %colormap for plotting

%calculate number of samples and probability correct for each
%condition
nc = accumarray([sIndex, cIndex, dxIndex], data.response);
ncc = accumarray([sIndex, cIndex, dxIndex], 1-data.response);
n = nc + ncc;
pc = nc ./ n; %this will be NaN

%Loop through each spacing and plot the psychometric functions

%Obtain predictions from the model
model_pred = nan(size(n));
for sNum = 1:length(sList)
    for cNum = 1:length(cList)
        if sum(n(sNum,cNum,:))
            model_pred(sNum, cNum, :) = ...
                MotionModel(p,struct(...
                    'spacing', sList(sNum)*ones(size(dxList)), ...
                    'content', cList(cNum)*ones(size(dxList)), ...
                    'dx', dxList));
        end
    end
end

%also fit cumulative normals per condition. They go into mu and
%sig, indexed by sNum and cNum.
norm_pred = nan(size(n));
for sNum = 1:length(sList)
    for cNum = 1:length(cList)
        %Initial parameters for the cumulative normal
        [~,pNorm.mu,pNorm.sig] = ...
            MotionModel(p,struct(...
                'spacing', sList(sNum), 'content', cList(cNum), 'dx', 0));
        pNorm.shutup = 'yes';

        %Find all trials for this s and c (and all dx)
        id = data.spacing == sList(sNum) & data.content == cList(cNum);
        if sum(id)
            %Fit the cumulative normal (with slope fixed by model fit?)
            pNormBest = fit(@fitCumNormal,pNorm,{'mu'}, ...
                            data.dx(id),data.response(id));
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
        h(cNum) = plot(dxList,squeeze(pc(sNum,cNum,:)),':',...
                       'Color',colList(cNum,:));

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
    [~,predMu,~] = ...
        MotionModel(p,struct('spacing', ...
                             sList, 'content', cList(cNum), 'dx', 0));
    plot(sList,predMu,'-','Color',colList(cNum,:));
end
xlabel('s');
ylabel('mu');
legend(h,num2str(cList));

subplot(2,1,2); %sig as a function of s and c
hold on

for cNum =1:length(cList)
    %Plot sig for the best fitting cumulative normal
    id = sig(:,cNum) >0;
    h(cNum) = plot(sList(id),sig(id,cNum),'o','Color',colList(cNum,:));
    %get and plot sig from the model fit
    [prob,predMu,predSig] = ...
        MotionModel(p,struct('spacing', sList, ...
                             'content', cList(cNum), 'dx',  0));
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

save modelResults.mat
if (exist('export', 'var') && export)
    for varname = who()'
        assignin('base', varname{1}, eval(varname{1}));
    end
end
