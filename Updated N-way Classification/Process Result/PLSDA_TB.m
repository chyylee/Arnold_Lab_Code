%% PLSDA.m
% To perform PLSDA. You can use Iterate_Pscripts.m to perform PLSDA on all
% workspaces (made using Make_xyblocks.m) with in a given folder.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Melissa Lemke, Arnold Lab, University of Michigan, Biomedical Engineering
% March 16th, 2018
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% For iteration 
% function PLSDA(path)
% load(path)

%% Use only 'For iteration' section OR 'For analysis...') NEVER both
%% For analysis one workspace at a time


function [cls_error,cls_error_CV] = PLSDA_TB(xblock,yblock,xnames,filettl,classes,yname,ttl,flag)


%% Modeling
% Eigenvector options for the plsda modeling function
opts = plsda('options');
opts.plots = 'none';
opts.discrim = 'yes';
opts.orthogonalize = 'off'; % If you turn on preVIP, you need to use OPLS-VIP
opts.preprocessing = 'Autoscale';
opts.confidencelimit = 0.9;

% The logic to determine which splits and blindsizes for crossvalidation
if length(xblock(:,1))<20
    splits = floor(length(xblock(:,1))/2);
else
    splits = 10;
end
blindsize = 1;

cvi = {'vet' splits blindsize};

% Eigenvector options for the crossval modeling function
cvopts = crossval('options');
cvopts.preprocessing = 2;
cvopts.discrim = 'yes';

% The initial model to decide how many LVs to use for VIP selection

plsda_model = plsda(xblock,yblock,length(xblock(1,:)),opts);
plsda_model = crossval(xblock,yblock,plsda_model,cvi,length(xblock(1,:)),cvopts);

cls_error = mean(plsda_model.classerrc)';
cls_error_CV = mean(plsda_model.classerrcv)';

% To plot the Cal and CV error of the initial model to choose components
plot(mean(plsda_model.classerrc)','-+')
hold on
plot(mean(plsda_model.classerrcv)','-+')
legend('Calibration Error','Crossvalidation Error')
xlabel('Latent Variable Number')
ylabel('Error Average')
hold off

if flag == 1

% Prompt used to choose the number of LVs for VIP selection
inputprompt={'Choose # of LVs: for VIP selection: '};
inputname='Choose Components';
numlines=1;
defaultans={'2'};
sortoptions=inputdlg(inputprompt,inputname,numlines,defaultans);
num_LVs = str2num(sortoptions{1});
close

% Prompt used to choose the number of VIP selection rounds
inputprompt={'Choose # of VIP selection rounds: '};
inputname='Choose # VIP';
numlines=1;
defaultans={'1'};
sortoptions=inputdlg(inputprompt,inputname,numlines,defaultans);
num_VIPs = str2num(sortoptions{1});


% The model to perform vip on

lowerror_model = plsda(xblock,yblock,num_LVs,opts);
lowerror_model = crossval(xblock,yblock,lowerror_model,cvi,num_LVs,cvopts);
close 

model_to_vip = lowerror_model;
if num_VIPs>0
for i = 1:num_VIPs
% VIP selection
vipscores = vip(model_to_vip);
index = find(sum(vipscores>=1,2));
vip_xblock = xblock(:,index);
vip_xnames = xnames(index);

% The VIP model used to determine the # of LVs in the final VIP model
hold off
lowerror_vip_model = plsda(vip_xblock,yblock,length(index),opts);
lowerror_vip_model = crossval(vip_xblock,yblock,lowerror_vip_model,cvi,length(index),cvopts);
model_to_vip = lowerror_vip_model;
close
% To plot the Cal and CV error of the lowerror_vip_model to choose components

hold off
plot(mean(lowerror_vip_model.classerrc)','-+')
hold on
plot(mean(lowerror_vip_model.classerrcv)','-+')
legend('Calibration Error','Crossvalidation Error')
xlabel('Latent Variable Number')
ylabel('Error Average')
hold off

% Prompt to choose number of LVs in the vip selected model
inputprompt={'Choose # of LVs in model: '};
inputname='Choose Components';
numlines=1;
defaultans={'2'};
sortoptions=inputdlg(inputprompt,inputname,numlines,defaultans);
num_LVs=str2double(sortoptions{1});

close
opts.plots='off';
cvopts.plots='off';

% You can decide whether to orthogonalize the final model based on model
% interpretability, but pre VIP selection ortho must be off
opts.orthogonalize = 'on';

% Final VIP model
vip_model = plsda(vip_xblock,yblock,num_LVs,opts);
vip_model = crossval(vip_xblock,yblock,vip_model,cvi,num_LVs,cvopts);

model = vip_model;
end
else
    model = lowerror_model;
    vip_xblock = xblock;
    vip_xnames = xnames;
    end
 
calerror = mean(model.classerrc);
cverror = mean(model.classerrcv);

% 'Fake' LV2 model, in case you want 2D scores plot with a 1 LV model
if num_LVs == 1
    fakeLV2model = plsda(vip_xblock,yblock,2,opts);
    fakeLV2model = crossval(vip_xblock,yblock,fakeLV2model,cvi,2,cvopts);
end

%% Figure 
% Prompt to choose the # of LVs in the scores and loadings plot
inputprompt={'Choose # of LVs in scores plot: ';'Choose # of LVs in loadings plot: '};
inputname='FOR FIGURE';
numlines=1;
defaultans={'2',char(string(num_LVs))};
sortoptions=inputdlg(inputprompt,inputname,numlines,defaultans);
scores_num_LVs = str2num(sortoptions{1});
loads_num_LVs = str2num(sortoptions{2});

% Logic determining which scores loadings and axis labels to use in the
% figure. This depends on which model was actually chosen and then which
% dimensions were chosen for the figure.
LV1_scores = model.loads{1}(:,1);
LV1_loads = model.loads{2}(:,1);
scores_xlabel = strcat('Scores on LV1 (',num2str(model.ssq(1,2),4),'%)');
loads_xlabel = strcat('Loadings on LV1 (',num2str(model.ssq(1,2),4),'%)');

if num_LVs > 1 % 2 LV model
    LV2_scores = model.loads{1}(:,2);
    LV2_loads = model.loads{2}(:,2);
    scores_ylabel = strcat('Scores on LV2 (',num2str(model.ssq(2,2),4),'%)');
    loads_ylabel = ['Loadings on LV2 (',num2str(model.ssq(2,2),4),'%)'];
else % 1 LV model
    LV2_scores = [];
    LV2_loads = [];
    scores_ylabel = [];
    loads_ylabel = [];
    if scores_num_LVs == 2 % if model is 1 LV but you want the 2 LV scores plot
        LV1_scores = fakeLV2model.loads{1}(:,1);
        LV2_scores = fakeLV2model.loads{1}(:,2);
        scores_xlabel = strcat('Scores on LV1 (',num2str(fakeLV2model.ssq(1,2),4),'%)');
        scores_ylabel = strcat('Scores on LV2 (',num2str(fakeLV2model.ssq(2,2),4),'%)');
    end
end

% The text displayed between the scores and loadings plot at the bottom
LV = 1:num_LVs;
error = [LV(num_LVs),calerror(num_LVs),cverror(num_LVs)];
errortext = sprintf('         Cal Error  Cv Error\n LV%1.0f   %1.4f     %1.4f\n',error');

% The title for the workspace and figure saved into the current directory
modelfilettl = strcat(datestr(today()),strtok(extractAfter(filettl,11),'.'),'_PLSDA_Model.mat');

% Calls the ScoresandLoadingsPlot function to make the figure
ScoresandLoadingsPlot(LV1_scores, LV2_scores, [],LV1_loads, LV2_loads,[],...
    yblock, scores_xlabel, scores_ylabel, [], loads_xlabel, loads_ylabel, [],...
    classes,  vip_xnames, "", ttl, string(modelfilettl), errortext,0,...
    scores_num_LVs, loads_num_LVs)
hold off

% Calls the ScoresandLoadingsPlot function if you want to have two figures
% with 1D loadings plots for both PC1 and PC2 (i.e. creates the PC2
% loadings figure)
if scores_num_LVs == 2 && loads_num_LVs == 1 && num_LVs == 2
    ScoresandLoadingsPlot(LV1_scores, LV2_scores, [], LV2_loads,[], [],...
        yblock, scores_xlabel, scores_ylabel, [], loads_ylabel, [], [],...
        classes,  vip_xnames, "", ttl, strcat(strtok(modelfilettl,'.'),'_LV2.mat'),...
        errortext,0, scores_num_LVs, loads_num_LVs);
    hold off
end

% Saves the entire workspace, including figures to the current directory
save(string(modelfilettl))
end
end