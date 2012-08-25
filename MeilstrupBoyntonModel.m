% AnalyzePeter.m

% Notes 
% s1 is hard to fit all without c0 = -.1
% s2 is OK
% s3 is hard.  No dependency of mu on c
% s4 is pretty good
% s5 is sparse and OK except for s=2.62
% s6 is unfittable
% s7 is good
% s8 is pretty good
% s9 is OK
% s10 (pbm) is pretty good
% s11 is sparse
% s12 is pretty good
% s13 is sparse but OK anywa


% Which data sets to fit
expTypes = {'content','spacing','all'};
expNum = 2;  %3 fits both types together

subjectNum  =10;  %s10 is pbm

%Initial parameters
switch subjectNum
    case 1
        p.sig0 = .065;
        p.sigk = .75;
        p.siga = 3;
        
        p.mua = 2;
        p.mukc = 1.5;
        p.muks = 1.5;
        p.mu0 = -.1;
        
    case 3
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
if ~exist('contentData','var')
    load data
end
% Copy the approprate structure to 'data'
eval(sprintf('data = %sData;',expTypes{expNum}));

%% Pull out the relevant trials and fit the model to the data
% List of all subjects
subjectList = unique(data.subject);
nSubjects = length(subjectList);

% Find trials for the current subject
id = strcmp(data.subject,subjectList{subjectNum});

% Pull out all conditions and responses for all trials for this subject
s = data.target_spacing(id);         
c = data.folded_direction_content(id);
dx = data.folded_displacement(id);
response = 1-data.folded_response_with_carrier(id);  %note the 1-response.

%Fit the model (uncomment the next line to see initial parameter predictions)
p = fit('fitMotionModel',p,freeParams,s,c,dx,response);
p

%% Plot each psychometeric function and model prediction
close all

%Find the list of parameters used in the experiments using 'unique'
eccentrities = unique(data.eccentricity(id));  %not used, always one eccentricity (6.667 deg)
sList = unique(s);   
dxList = unique(dx);
cList = unique(c);

colList = cool(length(cList));  %colormap for plotting

%Loop through each spacing and plot the psychometric functions
clear n pc mu sig
for sNum = 1:length(sList)
    clear h
    figure(sNum)
    clf
    hold on
    title(sprintf('subject %s, s = %5.2f',subjectList{subjectNum},sList(sNum)));
    % loop through the c's, generating separate psychometric functions
    for cNum = 1:length(cList)
        % loop through the dx's
        for dxNum = 1:length(dxList)
            %Find all trials for this set of s, c and dx
            id = s==sList(sNum) & c == cList(cNum) & dx == dxList(dxNum);
            n(sNum,cNum,dxNum) = sum(id);  %number of trials
            pc(sNum,cNum,dxNum) = mean(response(id)); %percent clockwise
            if n(sNum,cNum,dxNum)  %plot that point, scaled in size if there are any trials here
                plot(dxList(dxNum),pc(sNum,cNum,dxNum),'o',...
                    'MarkerSize',n(sNum,cNum,dxNum)/20+5,'MarkerFaceColor',colList(cNum,:));
            end
        end
        %plot a dashed line through the points
        h(cNum) = plot(dxList,squeeze(pc(sNum,cNum,:)),':','Color',colList(cNum,:));
        
        %plot the model prediction as a solid line
        
        if sum(n(sNum,cNum,:))
            pred = MotionModel(p,sList(sNum)*ones(size(dxList)),cList(cNum)*ones(size(dxList)),dxList);
            plot(dxList,pred,'-','Color',colList(cNum,:));
        end
        
        %fit the psychometric functions separately
        
        %Initial parameters for the cumulative normal
        
            [prob,predMu,predSig] = MotionModel(p,sList(sNum),cList(cNum),0);

        
        pNorm.sig =  predSig;       
        pNorm.mu = predMu;
        pNorm.shutup = 'yes'; 
        
        %Find all trials for this s and c (and all dx)
        id = s==sList(sNum) & c == cList(cNum);
        if sum(id)
            %Fit the cumulative normal 
            pNormBest = fit('fitCumNormal',pNorm,{'mu'},dx(id),response(id));
            %Save the best-fitting parameters
            mu(sNum,cNum) = pNormBest.mu;
            sig(sNum,cNum) = pNormBest.sig;
            
            %Plot the best fitting cumulative normal
            pred = normcdf(dxList,mu(sNum,cNum),sig(sNum,cNum));
            %plot(dxList,pred,'-.','Color',colList(cNum,:));
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
switch(expNum)
    case 1
        tile(2,2)
    case {2, 3}
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





